#include "../analysis_tools.C"

void plot_ramp(TString folder, TString output_file="", TString current_filename=""){

	TGraphErrors* g_ramp = new TGraphErrors();
	TGraph* g_rampA = new TGraph();
	TGraph* g_rampB = new TGraph();
	TGraph* g_rampAB = new TGraph();
	TGraph* g_ramp_norm = new TGraph();

	TGraph* g_diff_vs_C = new TGraph();

	map<int,double> map_time_current;
	bool know_current_info = false;

	if (current_filename != ""){
		know_current_info = true;
		ifstream file_current;
		file_current.open(current_filename);
		if (!file_current.is_open()){
			cout<<"cannot open "<<current_filename<<"\n";
			return;
		}

		char buff[256];
		file_current.getline(buff,256); //Header
		double time, current;
		char comma;
		while (!file_current.eof()){

			file_current.getline(buff,256);
			stringstream iss(buff);
			string substr;
			getline(iss,substr,',');
			if (substr == "") continue;

			int yy, mm, dd, hh, mi, ss;
    		if (!(sscanf(substr.c_str(), "%d-%d-%dT%d:%d:%d", &yy, &mm, &dd, &hh, &mi, &ss) == 6)){
        		cout<<"Can't extract time from string "<<substr<<"!\n";
        		throw std::runtime_error{"failed to parse time string"};
    		}
    		TDatime this_time = TDatime(yy,mm,dd,hh,mi,ss);
    		int time_stamp = this_time.Convert();

			getline(iss,substr,',');
			if (substr == "") continue;
			current = stof(substr);

			if (file_current.eof()) break;
			map_time_current[time_stamp] = current;
		}
		file_current.close();
	}

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	double ABref = 0;
	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%100==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_A = traces[map_varnames["Channel A"]];
		vector<double> trace_B = traces[map_varnames["Channel B"]];
		vector<double> trace_diff = traces[map_varnames["average(B-A)"]];
		vector<double> trace_C = traces[map_varnames["Channel C"]];

		TDatime datetime = getFileTime(fname);
		int time_stamp = datetime.Convert();

		double current = fi;
		if (current_filename != ""){
			current = map_time_current[time_stamp];
		}
		
		double A_avg = 0;
		double B_avg = 0;
		double C_avg = 0;
		double ABsum = 0;
		double ABdiff_avg = 0;
		double ABdiff_avg_squared = 0;
		double Npoints = trace_time.size();
		for (int i=0; i<Npoints; i++){
			A_avg += trace_A[i];
			B_avg += trace_B[i];
			C_avg += trace_C[i];
			ABdiff_avg += trace_diff[i];
			ABdiff_avg_squared += trace_diff[i]*trace_diff[i];
		}
		A_avg /= Npoints;
		B_avg /= Npoints;
		C_avg /= Npoints;
		ABdiff_avg /= Npoints;
		ABdiff_avg_squared /= Npoints;
		double ABerr = sqrt(ABdiff_avg_squared - ABdiff_avg*ABdiff_avg);
		ABsum = A_avg + B_avg;
		if (fi==0) ABref = ABsum;

		int ipoint = g_ramp->GetN();
		g_ramp->SetPoint(ipoint,current,ABdiff_avg);
		g_ramp->SetPointError(ipoint,0,ABerr);
		g_rampA->SetPoint(ipoint,current,A_avg);
		g_rampB->SetPoint(ipoint,current,B_avg);
		g_rampAB->SetPoint(ipoint,current,A_avg+B_avg);
		g_ramp_norm->SetPoint(ipoint,current,ABdiff_avg*ABref/ABsum);

		if (abs(C_avg) < 1000){
			g_diff_vs_C->SetPoint(g_diff_vs_C->GetN(),B_avg-A_avg,C_avg);
		}
	}

	new TCanvas();
	g_ramp->SetName("Ramp");
	g_ramp->SetTitle("Ramp");
	g_ramp->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_ramp->GetYaxis()->SetTitle("B-A [V]");
	g_ramp->SetMarkerStyle(20);
	g_ramp->Draw("APL");
	gPad->SetGridy();

	new TCanvas();
	g_rampA->SetName("A");
	g_rampA->SetTitle("A");
	g_rampB->SetName("B");
	g_rampB->SetTitle("B");
	g_rampAB->SetName("AB");
	g_rampAB->SetTitle("AB");
	g_rampA->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampA->GetYaxis()->SetTitle("A [V]");
	g_rampB->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampB->GetYaxis()->SetTitle("B [V]");
	g_rampAB->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampAB->GetYaxis()->SetTitle("A+B [V]");
	g_rampA->GetYaxis()->SetRangeUser(0,15);
	g_rampB->GetYaxis()->SetRangeUser(0,15);
	g_rampAB->GetYaxis()->SetRangeUser(0,15);
	g_rampA->SetMarkerStyle(20);
	g_rampB->SetMarkerStyle(20);
	g_rampAB->SetMarkerStyle(20);
	g_rampA->SetMarkerColor(kBlue);
	g_rampB->SetMarkerColor(kRed);
	g_rampAB->SetMarkerColor(kBlack);
	g_rampA->SetLineColor(kBlue);
	g_rampB->SetLineColor(kRed);
	g_rampAB->SetLineColor(kBlack);
	g_rampA->Draw("APL");
	g_rampB->Draw("PL");
	g_rampAB->Draw("PL");
	gPad->SetGridy();


	new TCanvas();
	g_ramp_norm->SetName("Ramp_norm");
	g_ramp_norm->SetTitle("Ramp normalized");
	g_ramp_norm->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_ramp_norm->GetYaxis()->SetTitle("(B-A)/(A+B)");
	g_ramp_norm->SetMarkerStyle(20);
	g_ramp_norm->Draw("APL");
	gPad->SetGridy();


	new TCanvas();
	g_diff_vs_C->SetMarkerStyle(20);
	g_diff_vs_C->Draw("AP");
	g_diff_vs_C->Fit("pol1");


	if (output_file != ""){
		cout<<"Creating "<<output_file<<"\n";
		TFile* fout = new TFile(output_file,"recreate");
		g_ramp->Write();
		g_rampA->Write();
		g_rampB->Write();
		g_rampAB->Write();
		g_ramp_norm->Write();
		fout->Write();
		fout->Close();
	}

}
#include "../analysis_tools.C"

void plot_rampup(TString folder, TString output_file, TString current_filename=""){

	TGraphErrors* g_rampup = new TGraphErrors();
	TGraph* g_rampupA = new TGraph();
	TGraph* g_rampupB = new TGraph();


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

	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%100==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_A = traces[map_varnames["Channel A"]];
		vector<double> trace_B = traces[map_varnames["Channel B"]];
		vector<double> trace_diff = traces[map_varnames["average(B-A)"]];

		TDatime datetime = getFileTime(fname);
		int time_stamp = datetime.Convert();

		double current = fi;
		if (current_filename != ""){
			current = map_time_current[time_stamp];
		}
		
		double A_avg = 0;
		double B_avg = 0;
		double ABdiff_avg = 0;
		double ABdiff_avg_squared = 0;
		double Npoints = trace_time.size();
		for (int i=0; i<Npoints; i++){
			A_avg += trace_A[i];
			B_avg += trace_B[i];
			ABdiff_avg += trace_diff[i];
			ABdiff_avg_squared += trace_diff[i]*trace_diff[i];
		}
		A_avg /= Npoints;
		B_avg /= Npoints;
		ABdiff_avg /= Npoints;
		ABdiff_avg_squared /= Npoints;
		double ABerr = sqrt(ABdiff_avg_squared - ABdiff_avg*ABdiff_avg);

		int ipoint = g_rampup->GetN();
		g_rampup->SetPoint(ipoint,current,ABdiff_avg);
		g_rampup->SetPointError(ipoint,0,ABerr);
		g_rampupA->SetPoint(ipoint,current,A_avg);
		g_rampupB->SetPoint(ipoint,current,B_avg);
	}

	new TCanvas();
	//g_rampup->Sort();
	g_rampup->SetName("Rampup");
	g_rampup->SetTitle("Ramp up");
	g_rampup->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampup->GetYaxis()->SetTitle("B-A [V]");
	g_rampup->SetMarkerStyle(20);
	g_rampup->Draw("APL");

	new TCanvas();
	//g_rampupA->Sort();
	//g_rampupB->Sort();
	g_rampupA->SetName("A");
	g_rampupA->SetTitle("A");
	g_rampupB->SetName("B");
	g_rampupB->SetTitle("B");
	g_rampupA->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampupA->GetYaxis()->SetTitle("B-A [V]");
	g_rampupB->GetXaxis()->SetTitle(know_current_info?"Current [A]":"File number");
	g_rampupB->GetYaxis()->SetTitle("B-A [V]");
	g_rampupA->SetMarkerStyle(20);
	g_rampupB->SetMarkerStyle(20);
	g_rampupA->SetMarkerColor(kBlue);
	g_rampupB->SetMarkerColor(kRed);
	g_rampupA->SetLineColor(kBlue);
	g_rampupB->SetLineColor(kRed);
	g_rampupA->Draw("APL");
	g_rampupB->Draw("PL");

	cout<<"Creating "<<output_file<<"\n";
	TFile* fout = new TFile(output_file,"recreate");
	g_rampup->Write();
	g_rampupA->Write();
	g_rampupB->Write();
	fout->Write();
	fout->Close();


}
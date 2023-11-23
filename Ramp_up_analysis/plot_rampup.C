#include <dirent.h>

void plot_rampup(TString folder, TString output_file, TString current_filename=""){

	TGraph* g_trace = new TGraph();
	TGraph* g_traceA = new TGraph();
	TGraph* g_traceB = new TGraph();
	TGraphErrors* g_rampup = new TGraphErrors();
	TGraphErrors* g_rampupA = new TGraphErrors();
	TGraphErrors* g_rampupB = new TGraphErrors();



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
			stringstream ss(buff);
			string substr;
			getline(ss,substr,',');
			if (substr == "") continue;
			std::tm t{};
			std::istringstream iss(substr);
			iss >> std::get_time(&t,"%Y-%m-%dT%H:%M:%S");
			if (iss.fail()) {
				cout<<substr<<"\n";
			    throw std::runtime_error{"failed to parse time string"};
			}   
			int time_stamp = mktime(&t);
			getline(ss,substr,',');
			if (substr == "") continue;
			current = stof(substr);

			if (file_current.eof()) break;
			map_time_current[time_stamp] = current;
		}
		file_current.close();
	}


	DIR* dir = opendir(folder.Data());
	struct dirent* dirfile;
	cout<<"Reading folder "<<folder<<"\n";
	
	vector<TString> files;
	while((dirfile = readdir(dir)) != NULL){
		if (dirfile->d_type != DT_REG) continue;
		TString fname = dirfile->d_name;
		if (!fname.EndsWith("csv")) continue;
		files.push_back(fname);
	}

	sort(files.begin(),files.end());

	int Nfiles = files.size();

	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%100==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		ifstream file;
		TString filename = Form("%s/%s",folder.Data(),fname.Data());
		file.open(filename.Data());
		if (!file.is_open()){
			cout<<"cannot open "<<filename<<"\n";
			return;
		}
		char buff[256];
		file.getline(buff,256); //Header, names

		stringstream ss(buff);
		int Nvar = 0;
		double Apos = -1;
		double Bpos = -1;
		double ABpos = -1;
		while(ss.good()){
			string substr;
			getline(ss,substr,',');
			if (substr.find("average(B-A)") != string::npos){
				ABpos = Nvar;
			} else if (substr.find("A") != string::npos){
				Apos = Nvar;
			} else if (substr.find("B") != string::npos){
				Bpos = Nvar;
			}
			Nvar++;
		}

		file.getline(buff,256); //Header, units
		file.getline(buff,256); //Header, blank row

		TString file_datetime = fname;
		file_datetime.ReplaceAll("FD_","");
		file_datetime.ReplaceAll("Ramp_","");
		file_datetime.ReplaceAll(".csv","");
		std::tm t{};
		std::istringstream iss(file_datetime.Data());
		iss >> std::get_time(&t, "%m_%d_%Y %H_%M_%S");
		if (iss.fail()) {
			cout<<file_datetime<<"\n";
		    throw std::runtime_error{"failed to parse time string"};
		}   
		int time_stamp = mktime(&t);

		double current = fi;

		if (current_filename != ""){
			current = map_time_current[time_stamp];
		}
		
		double time, A, B, C, avgAB;
		char comma;
		g_trace->Set(0);
		g_traceA->Set(0);
		g_traceB->Set(0);
		while (!file.eof()){
			double var=0;
			for (int i=0; i<Nvar; i++){
				if (i>0) file>>comma;
				file>>var;

				if (i==0){time = var;}
				if (i==Apos){A = var;}
				if (i==Bpos){B = var;}
				if (i==ABpos){avgAB = var;}
			}

			if (file.eof()) break;

			g_trace->AddPoint(time,avgAB);
			g_traceA->AddPoint(time,A);
			g_traceB->AddPoint(time,B);
		}
		file.close();


		TFitResultPtr fit_flat = g_trace->Fit("pol0","QNS");
		TFitResultPtr fit_flatA = g_traceA->Fit("pol0","QNS");
		TFitResultPtr fit_flatB = g_traceB->Fit("pol0","QNS");

		double ABval = fit_flat->Parameter(0);
		double ABerr = fit_flat->ParError(0);
		double Aval = fit_flatA->Parameter(0);
		double Aerr = fit_flatA->ParError(0);
		double Bval = fit_flatB->Parameter(0);
		double Berr = fit_flatB->ParError(0);


		g_rampup->AddPoint(current,ABval);
		g_rampupA->AddPoint(current,Aval);
		g_rampupB->AddPoint(current,Bval);
		int ipoint = g_rampup->GetN()-1;
		g_rampup->SetPointError(ipoint,0,ABerr);
		g_rampupA->SetPointError(ipoint,0,Aerr);
		g_rampupB->SetPointError(ipoint,0,Berr);

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
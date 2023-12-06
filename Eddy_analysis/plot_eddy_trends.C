#include <dirent.h>
#include "../analysis_tools.C"


void plot_eddy_trends(TString folder, TString filename_out){


	TFile* fout = new TFile(filename_out,"recreate");

	TProfile* p_trace = new TProfile("p_trace","",195315,0,100);
	TGraph* g_trace = new TGraph();

	TGraph* g_trend_ABdiff = new TGraph();
	TGraph* g_trend_Cstddev = new TGraph();
	TGraph* g_trend_avgCstddev = new TGraph();

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
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		vector<pair<double,double>> trace;

		ifstream file;
		TString filename = Form("%s/%s",folder.Data(),fname.Data());
		file.open(filename.Data());
		if (!file.is_open()){
			cout<<"cannot open "<<filename<<"\n";
			return;
		}


		TString file_datetime = fname;
		file_datetime.ReplaceAll("SD_","");
		file_datetime.ReplaceAll(".csv","");
		std::tm tm{};
		std::istringstream iss(file_datetime.Data());
		iss >> std::get_time(&tm, "%m_%d_%Y %H_%M_%S");
		if (iss.fail()) {
			cout<<file_datetime<<"\n";
		    throw std::runtime_error{"failed to parse time string"};
		}

		TDatime datetime = TDatime(tm.tm_year+1900, tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

		char buff[256];
		file.getline(buff,256); //Header, names
		stringstream ss(buff);
		int Nvar = 0;
		double Apos = -1;
		double Bpos = -1;
		double Cpos = -1;
		double avgCpos = -1;
		while(ss.good()){
			string substr;
			getline(ss,substr,',');
			if (substr.find("average(C)") != string::npos){
				avgCpos = Nvar;
			} else if (substr.find("A") != string::npos){
				Apos = Nvar;
			} else if (substr.find("B") != string::npos){
				Bpos = Nvar;
			} else if (substr.find("C") != string::npos){
				Cpos = Nvar;
			}
			Nvar++;
		}

		file.getline(buff,256); //Header, units
		file.getline(buff,256); //Header, blank row

		double time, A, B, C, avgC;
		double avgA = 0;
		double avgB = 0;
		double avgC_ = 0;
		double avgC2_ = 0;
		double avgavgC_ = 0;
		double avgavgC2_ = 0;
		double Cstddev = 0;
		double avgCstddev = 0;
		char comma;
		double Nlines = 0;
		double Nlines_4ms = 0;
		while (!file.eof()){
			double var=0;
			for (int i=0; i<Nvar; i++){
				if (i>0) file>>comma;
				file>>var;

				if (i==0){time = var;}
				if (i==Apos){A = var;}
				if (i==Bpos){B = var;}
				if (i==Cpos){C = var;}
				if (i==avgCpos){avgC = var;}
			}

			if (file.eof()) break;

			avgA += A;
			avgB += B;
			avgC_ += C;
			avgC2_ += C*C;

			if (time < 4){
				avgavgC_ += avgC;
				avgavgC2_ += avgC*avgC;
				Nlines_4ms++;
			}

			trace.push_back(make_pair(time,avgC));
			Nlines++;
		}
		file.close();

		avgA /= Nlines;
		avgB /= Nlines;
		avgC_ /= Nlines;
		avgC2_ /= Nlines;
		Cstddev = sqrt(avgC2_ - avgC_*avgC_);
		avgavgC_ /= Nlines;
		avgavgC2_ /= Nlines;
		avgCstddev = sqrt(avgavgC2_ - avgavgC_*avgavgC_);
		
		g_trend_ABdiff->SetPoint(g_trend_ABdiff->GetN(),datetime.Convert(),avgB-avgA);
		g_trend_Cstddev->SetPoint(g_trend_Cstddev->GetN(),datetime.Convert(),Cstddev);
		g_trend_avgCstddev->SetPoint(g_trend_avgCstddev->GetN(),datetime.Convert(),avgCstddev);


		g_trace->Set(0);
		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_trace->SetPoint(g_trace->GetN(),time,avgC);
		}

		TFitResultPtr res = g_trace->Fit("pol0","QNS","",2.0,12.0);
		if (res>=0){
			double baseline = res->Parameter(0);
			for (auto p : trace){
				time = p.first;
				avgC = p.second;
				p_trace->Fill(time,avgC-baseline);
			}
		}
	}


	g_trend_ABdiff->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_ABdiff->GetXaxis()->SetTimeOffset(0);
	g_trend_ABdiff->GetXaxis()->SetTimeDisplay(1);
	g_trend_Cstddev->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_Cstddev->GetXaxis()->SetTimeOffset(0);
	g_trend_Cstddev->GetXaxis()->SetTimeDisplay(1);
	g_trend_avgCstddev->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_avgCstddev->GetXaxis()->SetTimeOffset(0);
	g_trend_avgCstddev->GetXaxis()->SetTimeDisplay(1);

	new TCanvas();
	p_trace->GetXaxis()->SetTitle("Time [ms]");
	p_trace->GetYaxis()->SetTitle("Trace [mV]");
	p_trace->SetMarkerStyle(8);
	p_trace->Draw("HIST");

	new TCanvas();
	g_trend_ABdiff->SetTitle("B-A;File;B-A [V]");
	g_trend_ABdiff->SetMarkerStyle(20);
	g_trend_ABdiff->Draw("APL");
	
	new TCanvas();
	g_trend_Cstddev->SetTitle("C stddev;File;C stddev [V]");
	g_trend_Cstddev->SetMarkerStyle(20);
	g_trend_Cstddev->Draw("APL");

	new TCanvas();
	g_trend_avgCstddev->SetTitle("avgC stddev before kick1;File;avgC stddev [mV]");
	g_trend_avgCstddev->SetMarkerStyle(20);
	g_trend_avgCstddev->Draw("APL");


	fout->Write();
	fout->Close();
}

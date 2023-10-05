#include <dirent.h>


void plot_eddy_simple(TString folder){


	TProfile* p_trace = new TProfile("p_trace","",195315,0,100);
	TGraph* g_trace = new TGraph();

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
		char buff[256];
		file.getline(buff,256); //Header, names

		stringstream ss(buff);
		int Nvar = 0;
		double Cpos = -1;
		while(ss.good()){
			string substr;
			getline(ss,substr,',');
			//cout<<Nvar<<" : "<<substr<<"\n";
			if (substr.find("average(C)") != string::npos){
				Cpos = Nvar;
			}
			Nvar++;
		}

		file.getline(buff,256); //Header, units
		file.getline(buff,256); //Header, blank row

		double time, A, B, C, avgC;
		char comma;
		while (!file.eof()){
			double var=0;
			for (int i=0; i<Nvar; i++){
				if (i>0) file>>comma;
				file>>var;

				if (i==0){time = var;}
				if (i==Cpos){avgC = var;}
			}

			if (file.eof()) break;

			trace.push_back(make_pair(time,avgC));
		}
		file.close();

		g_trace->Set(0);
		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_trace->AddPoint(time,avgC);
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

	new TCanvas();
	p_trace->GetXaxis()->SetTitle("Time [ms]");
	p_trace->GetYaxis()->SetTitle("Trace [mV]");
	p_trace->SetMarkerStyle(8);
	p_trace->Draw("HIST");
}

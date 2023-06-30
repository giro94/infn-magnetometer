#include <dirent.h>


void plot_HWPscan(TString folder, double hwp_start=0, double hwp_step=5){


	TGraph** g_traces = new TGraph* [100];
	TGraph* g_scan = new TGraph();

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
		if (fi%10==0) cout<<"Reading file "<<fname<<"\n";

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

			//file>>time>>comma>>A>>comma>>B>>comma>>C>>comma>>avgC;
			//file>>time>>comma>>C>>comma>>avgC;
			trace.push_back(make_pair(time,avgC));
		}
		file.close();

		g_traces[fi] = new TGraph();
		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_traces[fi]->AddPoint(time,avgC);
		}

		double miny = 1e6;
		double maxy = -1e6;
		for (auto p : trace){
			time = p.first;
			avgC = p.second;

			if (avgC < miny) miny = avgC;
			if (avgC > maxy) maxy = avgC;
		}

		double HWPangle = hwp_start + fi*hwp_step;
		double peak = maxy-miny;
		g_scan->AddPoint(HWPangle,peak);
	}

	new TCanvas();
	g_scan->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_scan->GetYaxis()->SetTitle("Signal amplitude [V]");
	g_scan->SetMarkerStyle(20);
	g_scan->Draw("APL");

	gStyle->SetPalette(kRainBow);
	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_traces[i]->GetXaxis()->SetTitle("Time [ms]");
		g_traces[i]->GetYaxis()->SetTitle("Trace [V]");
		g_traces[i]->SetLineWidth(2);
		g_traces[i]->Draw(i==0?"AL PLC":"L PLC");
	}
}
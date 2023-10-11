#include <dirent.h>

void plot_HWPscan(TString folder, vector<double> hwp_angles){

	bool do_sort = true;

	int Nangles = hwp_angles.size();

	TGraph** g_traces = new TGraph* [Nangles];
	TGraph** g_baselines = new TGraph* [Nangles];
	TGraph** g_traces_norm = new TGraph* [Nangles];
	TGraphErrors* g_scan = new TGraphErrors();
	TGraph* g_stddev = new TGraph();
	TGraph* g_snr = new TGraph();

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

	if (Nfiles != Nangles){
		cout<<"Errors! Looking for "<<Nangles<<" HWP angles, but found "<<Nfiles<<" files!\n";
		return;
	}


	double baseline_fit_start = 0.1;
	double baseline_fit_end = 0.5;
	double blumlein_fit_start = 0.65;
	double blumlein_fit_end = 0.90; 

	TF1* f_baseline = new TF1("f_baseline","[0]",0,2);
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",0,2);

	for (int fi=0; fi<Nangles; fi++){

		TString fname = files[fi];
		cout<<"Reading file "<<fname<<" (angle "<<hwp_angles[fi]<<")\n";

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

			trace.push_back(make_pair(time,avgC));
		}
		file.close();

		g_traces[fi] = new TGraph();
		g_baselines[fi] = new TGraph();
		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_traces[fi]->AddPoint(time,avgC);

			if (time >= baseline_fit_start && time < baseline_fit_end){
				g_baselines[fi]->AddPoint(time,avgC);
			}
		}

		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baseline = g_traces[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);

		f_blumlein->SetParameters(60.0,0.5*(blumlein_fit_start+blumlein_fit_end),-1000.0);
		TFitResultPtr fit_blumlein = g_traces[fi]->Fit("f_blumlein","QS+","",blumlein_fit_start,blumlein_fit_end);

		double HWPangle = hwp_angles[fi];
		double baseline = fit_baseline->Parameter(0);
		double peak = fit_blumlein->Parameter(0);
		double baseline_err = fit_baseline->ParError(0);
		double peak_err = fit_blumlein->ParError(0);
		double amplitude = abs(peak-baseline);
		double baseline_stddev = g_baselines[fi]->GetRMS(2);
		double SNR = amplitude/baseline_stddev;

		g_stddev->AddPoint(HWPangle,baseline_stddev);
		g_snr->AddPoint(HWPangle,SNR);
		g_scan->AddPoint(HWPangle,amplitude);
		int ipoint = g_scan->GetN()-1;
		g_scan->SetPointError(ipoint,0,peak_err+baseline_err);


		g_traces_norm[fi] = new TGraph();
		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_traces_norm[fi]->AddPoint(time,avgC-baseline);
		}
	}

	new TCanvas();
	if (do_sort) g_snr->Sort();
	g_snr->SetName("HWPscan_SNR");
	g_snr->SetTitle("HWP scan SNR");
	g_snr->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_snr->GetYaxis()->SetTitle("Signal/Noise");
	g_snr->SetMarkerStyle(20);
	g_snr->SetLineStyle(kDashed);
	g_snr->Draw("APL");
	gPad->SetGridx();

	new TCanvas();
	if (do_sort) g_stddev->Sort();
	g_stddev->SetName("HWPscan_stddev");
	g_stddev->SetTitle("HWP scan StdDev");
	g_stddev->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_stddev->GetYaxis()->SetTitle("Baseline StdDev [mV]");
	g_stddev->SetMarkerStyle(20);
	g_stddev->Draw("APL");
	gPad->SetGridx();

	new TCanvas();
	if (do_sort) g_scan->Sort();
	g_scan->SetName("HWPscan");
	g_scan->SetTitle("HWP scan");
	g_scan->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_scan->GetYaxis()->SetTitle("Blumlein amplitude [mV]");
	g_scan->SetMarkerStyle(20);
	g_scan->Draw("APL");
	gPad->SetGridx();

	gStyle->SetPalette(kRainBow);
	new TCanvas();
	for (int i=0; i<Nangles; i++){
		g_traces[i]->SetName(Form("trace_%.1f",hwp_angles[i]));
		g_traces[i]->SetTitle(Form("Trace (hwp = %.1f)",hwp_angles[i]));
		g_traces[i]->GetYaxis()->SetRangeUser(-70,70);
		g_traces[i]->GetXaxis()->SetTitle("Time [ms]");
		g_traces[i]->GetYaxis()->SetTitle("Trace [mV]");
		g_traces[i]->Draw(i==0?"AL PLC":"L PLC");
	}

	new TCanvas();
	for (int i=0; i<Nangles; i++){
		g_traces_norm[i]->SetName(Form("trace_norm_%.1f",hwp_angles[i]));
		g_traces_norm[i]->SetTitle(Form("Normalized trace (hwp = %.1f)",hwp_angles[i]));
		g_traces_norm[i]->GetYaxis()->SetRangeUser(-70,70);
		g_traces_norm[i]->GetXaxis()->SetTitle("Time [ms]");
		g_traces_norm[i]->GetYaxis()->SetTitle("Trace [mV]");
		g_traces_norm[i]->Draw(i==0?"AL PLC":"L PLC");
	}


	TString outname = Form("%s_output.root",folder.Data());
	outname.ReplaceAll("../","");
	outname.ReplaceAll("/","__");
	outname.ReplaceAll(" ","_");
	cout<<"Creating "<<outname<<"\n";
	TFile* fout = new TFile(outname,"recreate");

	for (int i=0; i<Nangles; i++){
		g_traces[i]->Write();
	}
	g_scan->Write();
	g_stddev->Write();
	g_snr->Write();
	fout->Write();
	fout->Close();


}





void plot_HWPscan(TString folder, int N, double hwp_start=0, double hwp_step=5){

	vector<double> hwp_angles;
	for (int i=0; i<N; i++){
		hwp_angles.push_back(hwp_start+i*hwp_step);
	}
	plot_HWPscan(folder, hwp_angles);
}

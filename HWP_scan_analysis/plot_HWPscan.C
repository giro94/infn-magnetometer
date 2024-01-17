#include "../analysis_tools.C"

void plot_HWPscan(TString folder, vector<double> hwp_angles){

	bool do_sort = false;

	int Nangles = hwp_angles.size();

	TGraph** g_traces = new TGraph* [Nangles];
	TGraph** g_baselines = new TGraph* [Nangles];
	TGraph** g_traces_norm = new TGraph* [Nangles];
	TGraphErrors* g_scan = new TGraphErrors();
	TGraph* g_stddev = new TGraph();
	TGraph* g_snr = new TGraph();

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	if (Nfiles != Nangles){
		cout<<"Errors! Looking for "<<Nangles<<" HWP angles, but found "<<Nfiles<<" files!\n";
		return;
	}


	double baseline_fit_start;
	double baseline_fit_end;
	double blumlein_fit_start;
	double blumlein_fit_end;

	TF1* f_baseline = new TF1("f_baseline","[0]");
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");

	for (int fi=0; fi<Nangles; fi++){

		TString fname = files[fi];
		cout<<"Reading file "<<fname<<" (angle "<<hwp_angles[fi]<<")\n";

		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_avgC = traces[map_varnames["average(C)"]];

		g_traces[fi] = new TGraph();
		g_baselines[fi] = new TGraph();
		g_traces_norm[fi] = new TGraph();
		for (int i=0; i<trace_time.size(); i++){
			g_traces[fi]->SetPoint(i,trace_time[i],trace_avgC[i]);
		}

		//Find kick
		double first_kick_guess = -1;
		double firstkick_max = 20.0;
		double kick_trigger = -200.0;
		for (int i=1; i<trace_time.size(); i++){
			if (trace_avgC[i-1] > kick_trigger && trace_avgC[i] < kick_trigger){
				first_kick_guess = trace_time[i];
				break;
			}
		}
		if (first_kick_guess < 0){
			cout<<"Warning! Could not find first kick in file "<<fname<<". Skipping...\n";
			continue;
		} else if (first_kick_guess > firstkick_max){
			cout<<"Warning! First kick found at "<<first_kick_guess<<" in file "<<fname<<". Skipping...\n";
			continue;
		}

		//cout<<"Found kick at "<<first_kick_guess<<" ms\n";

		baseline_fit_start = first_kick_guess - 0.9;
		baseline_fit_end = first_kick_guess - 0.6;
		blumlein_fit_start = first_kick_guess - 0.43;
		blumlein_fit_end = first_kick_guess - 0.18;

		for (int i=0; i<trace_time.size(); i++){
			if (trace_time[i] >= baseline_fit_start && trace_time[i] < baseline_fit_end){
				g_baselines[fi]->SetPoint(i,trace_time[i],trace_avgC[i]);
			}
		}


		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baseline = g_traces[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);

		f_blumlein->SetParameters(60.0,0.5*(blumlein_fit_start+blumlein_fit_end),-1000.0);
		TFitResultPtr fit_blumlein = g_traces[fi]->Fit("f_blumlein","QS+","",blumlein_fit_start,blumlein_fit_end);

		if (fit_baseline<0 || fit_blumlein<0){
			cout<<"Bad fit!\n";
			continue;
		}

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


		for (int i=0; i<trace_time.size(); i++){
			g_traces_norm[fi]->SetPoint(i,trace_time[i],trace_avgC[i]-baseline);
		}
	}

	TCanvas* can_scan = new TCanvas("can_scan","",1800,600);
	can_scan->Divide(3,1);
	can_scan->cd(1);
	if (do_sort) g_scan->Sort();
	g_scan->SetName("HWPscan");
	g_scan->SetTitle("HWP scan");
	g_scan->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_scan->GetYaxis()->SetTitle("Blumlein amplitude [mV]");
	g_scan->SetMarkerStyle(20);
	g_scan->Draw("APL");
	gPad->SetGridx();

	can_scan->cd(2);
	if (do_sort) g_stddev->Sort();
	g_stddev->SetName("HWPscan_stddev");
	g_stddev->SetTitle("HWP scan StdDev");
	g_stddev->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_stddev->GetYaxis()->SetTitle("Baseline StdDev [mV]");
	g_stddev->SetMarkerStyle(20);
	g_stddev->Draw("APL");
	gPad->SetGridx();

	can_scan->cd(3);
	if (do_sort) g_snr->Sort();
	g_snr->SetName("HWPscan_SNR");
	g_snr->SetTitle("HWP scan SNR");
	g_snr->GetXaxis()->SetTitle("HWP angle [#circ]");
	g_snr->GetYaxis()->SetTitle("Signal/Noise");
	g_snr->SetMarkerStyle(20);
	g_snr->SetLineStyle(kDashed);
	g_snr->Draw("APL");
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

	TString outname = folder;
	outname.Remove(TString::kTrailing,'/');
    outname.Remove(0,outname.Last('/')+1);
    outname += "_output.root";
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


//Equally-spaced list of angles
void plot_HWPscan(TString folder, int N, double hwp_start=0, double hwp_step=5){

	vector<double> hwp_angles;
	for (int i=0; i<N; i++){
		hwp_angles.push_back(hwp_start+i*hwp_step);
	}
	plot_HWPscan(folder, hwp_angles);
}

//Angles from angles.txt file
void plot_HWPscan(TString folder){

	TString filepath = Form("%s/angles.txt",folder.Data());
	ifstream file_angles;
	file_angles.open(filepath.Data());
	if (!file_angles.is_open()){
		cout<<"cannot open "<<filepath<<"\n";
		throw std::runtime_error{"Cannot open angles.txt file, please specify angles as argument."};
	}

	vector<double> hwp_angles;
	while (!file_angles.eof()){
		double angle=0;
		file_angles>>angle;
		if (file_angles.eof()) break;
		hwp_angles.push_back(angle);
	}
	file_angles.close();

	plot_HWPscan(folder, hwp_angles);
}

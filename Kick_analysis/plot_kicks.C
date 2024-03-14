#include "../analysis_tools.C"

void plot_kicks(TString folder, TString output_file=""){

	TDatime oct8_k3ramp = TDatime(2023,10,8,10,38,0);
	TDatime oct9_k3ramp = TDatime(2023,10,9,15,25,0);

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	bool saveOutput=false;
	TFile* fout;
	if (output_file != ""){
		saveOutput=true;
		fout = new TFile(output_file,"recreate");
	}

	TGraph** g_kicks = new TGraph* [Nfiles];
	TGraph** g_kicksA = new TGraph* [Nfiles];
	TGraph** g_kicksB = new TGraph* [Nfiles];

	TGraph** g_kicks_linearized = new TGraph* [Nfiles];

	TGraph** g_kicks_aligned = new TGraph* [Nfiles];
	TGraph** g_kicks_normalized = new TGraph* [Nfiles];
	TGraph** g_kicks_diff = new TGraph* [Nfiles];

	double baseline_fit_start = 0.0;
	double baseline_fit_end = 1.0;
	double peak_fit_width = 0.05;

	TF1* f_baseline = new TF1("f_baseline","[0]");
	TF1* f_peak = new TF1("f_peak","[0]+[2]*(x-[1])*(x-[1])");

	int Nrepetitions = Nfiles/8;
	TGraphErrors** g_peaks = new TGraphErrors* [Nrepetitions];
	for (int i=0; i<Nrepetitions; i++){
		g_peaks[i] = new TGraphErrors();
	}
	TGraphErrors* g_peaks_time = new TGraphErrors();
	TGraphErrors* g_peaks_time_normalized = new TGraphErrors();
	TGraphErrors* g_A_time = new TGraphErrors();
	TGraphErrors* g_B_time = new TGraphErrors();
	TGraphErrors* g_AB_time = new TGraphErrors();

	double I = 0.8/sin(0.4);


	for (int fi=0; fi<Nfiles; fi++){
		
        TString fname = files[fi];
		TDatime datetime = getFileTime(fname);
		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		//Subtract kick ramp time
		if (fname.Index("10_08_2023") != -1){
			datetime = TDatime(datetime.Convert()-oct8_k3ramp.Convert());
		}
		if (fname.Index("10_09_2023") != -1){
			datetime = TDatime(datetime.Convert()-oct9_k3ramp.Convert());
		}

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_A = traces[map_varnames["Channel A"]];
		vector<double> trace_B = traces[map_varnames["Channel B"]];
		vector<double> trace_ABdiff = traces[map_varnames["average(B-A)"]];

		//Make sure we're looking at V and not mV
		vector<TString> headers = getHeaderUnits(filepath);
		if (headers[map_varnames["average(B-A)"]].Contains("mV")){
			for (int j=0; j<trace_ABdiff.size(); j++){
				trace_ABdiff[j] *= 0.001;
			}
		}

		g_kicks[fi] = new TGraph();
		g_kicksA[fi] = new TGraph();
		g_kicksB[fi] = new TGraph();
		g_kicks_linearized[fi] = new TGraph();
		for (int i=0; i<trace_time.size(); i++){
			g_kicks[fi]->SetPoint(i, trace_time[i], trace_ABdiff[i]);
			g_kicksA[fi]->SetPoint(i, trace_time[i], trace_A[i]);
			g_kicksB[fi]->SetPoint(i, trace_time[i], trace_B[i]);
			g_kicks_linearized[fi]->SetPoint(i, trace_time[i], asin(trace_ABdiff[i]/I));
		}


		int N = g_kicks[fi]->GetN();
		double peak_min = 0;
		double peak_max = 0;
		double peak_min_x = 0;
		double peak_max_x = 0;
		for (int j=0; j<N; j++){
			if (g_kicks[fi]->GetPointY(j) > peak_max){
				peak_max = g_kicks[fi]->GetPointY(j);
				peak_max_x = g_kicks[fi]->GetPointX(j);
			}
			if (g_kicks[fi]->GetPointY(j) < peak_min){
				peak_min = g_kicks[fi]->GetPointY(j);
				peak_min_x = g_kicks[fi]->GetPointX(j);
			}
		}

		if (abs(peak_min) > abs(peak_max)){ //invert trace
			for (int i=0; i<g_kicks[fi]->GetN(); i++){
				g_kicks[fi]->SetPointY(i,-g_kicks[fi]->GetPointY(i));
				g_kicksA[fi]->SetPointY(i,-g_kicksA[fi]->GetPointY(i));
				g_kicksB[fi]->SetPointY(i,-g_kicksB[fi]->GetPointY(i));
				g_kicks_linearized[fi]->SetPointY(i,-g_kicks_linearized[fi]->GetPointY(i));
			}
			peak_max_x = peak_min_x;
		}

		cout<<"Peak estimate: x="<<peak_max_x<<", y="<<peak_max<<"\n";

		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baseline = g_kicks[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);

		f_peak->SetParameters(0.8,peak_max_x,-1000.0);
		TFitResultPtr fit_peak = g_kicks[fi]->Fit("f_peak","QS+","",peak_max_x-peak_fit_width,peak_max_x+peak_fit_width);

		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baselineA = g_kicksA[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);

		f_peak->SetParameters(-0.8,peak_max_x,1000.0);
		TFitResultPtr fit_peakA = g_kicksA[fi]->Fit("f_peak","QS+","",peak_max_x-peak_fit_width,peak_max_x+peak_fit_width);

		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baselineB = g_kicksB[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);

		f_peak->SetParameters(0.8,peak_max_x,-1000.0);
		TFitResultPtr fit_peakB = g_kicksB[fi]->Fit("f_peak","QS+","",peak_max_x-peak_fit_width,peak_max_x+peak_fit_width);

		double baseline = fit_baseline->Parameter(0);
		double peak = fit_peak->Parameter(0);
		double baseline_err = fit_baseline->ParError(0);
		double peak_err = fit_peak->ParError(0);
		double amplitude = abs(peak-baseline);
		double amplitude_err = sqrt(peak_err*peak_err+baseline_err*baseline_err);
		double peak_time = fit_peak->Parameter(1);

		double baselineA = fit_baselineA->Parameter(0);
		double peakA = fit_peakA->Parameter(0);
		double baselineA_err = fit_baselineA->ParError(0);
		double peakA_err = fit_peakA->ParError(0);
		double amplitudeA = abs(peakA-baselineA);
		double amplitudeA_err = sqrt(peakA_err*peakA_err+baselineA_err*baselineA_err);

		double baselineB = fit_baselineB->Parameter(0);
		double peakB = fit_peakB->Parameter(0);
		double baselineB_err = fit_baselineB->ParError(0);
		double peakB_err = fit_peakB->ParError(0);
		double amplitudeB = abs(peakB-baselineB);
		double amplitudeB_err = sqrt(peakB_err*peakB_err+baselineB_err*baselineB_err);

		g_kicks_aligned[fi] = new TGraph();
		for (int j=0; j<N; j++){
			double x = g_kicks[fi]->GetPointX(j) - peak_time;
			double y = g_kicks[fi]->GetPointY(j) - baseline;
			g_kicks_aligned[fi]->SetPoint(j,x,y);
		}

		g_kicks_normalized[fi] = new TGraph();
		for (int j=0; j<N; j++){
			double x = g_kicks[fi]->GetPointX(j) - peak_time;
			double y = (g_kicks[fi]->GetPointY(j) - baseline) / amplitude; 
			g_kicks_normalized[fi]->SetPoint(j,x,y);
		}

		g_kicks_diff[fi] = new TGraph();
		for (int j=0; j<N; j++){
			double x = g_kicks_normalized[fi]->GetPointX(j);
			double y = g_kicks_normalized[fi]->GetPointY(j);
			double yref = 0;
			if (x >= g_kicks_normalized[0]->GetPointX(0) && x <= g_kicks_normalized[0]->GetPointX(N-1)){
				yref = g_kicks_normalized[0]->Eval(x); 
			}

			g_kicks_diff[fi]->SetPoint(j,x,y-yref);
		}

		int repetition_idx = (fi-fi%8)/8;
		int kick_idx = fi%8;
		g_peaks[repetition_idx]->SetPoint(kick_idx,kick_idx+1,amplitude);
		g_peaks[repetition_idx]->SetPointError(kick_idx,0,amplitude_err);

		g_peaks_time->SetPoint(fi,datetime.Convert(),amplitude);
		g_peaks_time->SetPointError(fi,0,amplitude_err);

		g_A_time->SetPoint(fi,datetime.Convert(),amplitudeA);
		g_A_time->SetPointError(fi,0,amplitudeA_err);
		g_B_time->SetPoint(fi,datetime.Convert(),amplitudeB);
		g_B_time->SetPointError(fi,0,amplitudeB_err);
		g_AB_time->SetPoint(fi,datetime.Convert(),amplitudeA+amplitudeB);
		g_peaks_time->SetPointError(fi,0,sqrt(amplitudeA_err*amplitudeA_err+amplitudeB_err*amplitudeB_err));

    }

    for (int i=0; i<g_peaks_time->GetN(); i++){
    	double norm_factor = 1.0/g_peaks_time->GetPointY(0);
    	g_peaks_time_normalized->SetPoint(i,g_peaks_time->GetPointX(i),g_peaks_time->GetPointY(i)*norm_factor);
    	g_peaks_time_normalized->SetPointError(i,g_peaks_time->GetErrorX(i),g_peaks_time->GetErrorY(i)*norm_factor);
    }


    TCanvas* can_kicks = new TCanvas("can_kicks","",1800,600);
    can_kicks->Divide(2,1);
    can_kicks->cd(1);
    for (int i=0; i<Nfiles; i++){
    	g_kicks[i]->SetTitle("All kicks");
    	g_kicks[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicks[i]->GetYaxis()->SetTitle("Trace [V]");
    	g_kicks[i]->GetXaxis()->SetRangeUser(0,6);
    	g_kicks[i]->Draw(i==0?"APL":"PL");
    }
    can_kicks->cd(2);
    for (int i=0; i<Nfiles; i++){
    	g_kicks_aligned[i]->SetTitle("Kicks aligned at t=0");
    	g_kicks_aligned[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicks_aligned[i]->GetYaxis()->SetTitle("Trace [V]");
    	g_kicks_aligned[i]->GetXaxis()->SetRangeUser(-1,5);
    	g_kicks_aligned[i]->Draw(i==0?"APL":"PL");
    }

    TCanvas* can_kicks_norm = new TCanvas("can_kicks_norm","",1800,600);
    can_kicks_norm->Divide(2,1);
    can_kicks_norm->cd(1);
    for (int i=0; i<Nfiles; i++){
    	g_kicks_normalized[i]->SetTitle("Normalized kicks");
    	g_kicks_normalized[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicks_normalized[i]->GetYaxis()->SetTitle("Trace [arb. u.]");
    	g_kicks_normalized[i]->GetXaxis()->SetRangeUser(-0.5,3.0);
    	g_kicks_normalized[i]->Draw(i==0?"APL":"PL");
    }
	gPad->SetGridy();
    can_kicks_norm->cd(2);
    for (int i=0; i<Nfiles; i++){
    	g_kicks_diff[i]->SetTitle("Normalized difference from first kick");
    	g_kicks_diff[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicks_diff[i]->GetYaxis()->SetTitle("Difference");
    	g_kicks_diff[i]->GetXaxis()->SetRangeUser(-0.5,3.0);
    	g_kicks_diff[i]->GetYaxis()->SetRangeUser(-0.05,0.05);
    	g_kicks_diff[i]->Draw(i==0?"APL":"PL");
    }
	gPad->SetGridy();


    TCanvas* can_amplitudes = new TCanvas("can_amplitudes","",1800,600);
    can_amplitudes->Divide(2,1);
    can_amplitudes->cd(1);
	g_peaks_time->SetMarkerStyle(20);
	g_peaks_time->SetTitle("Kick amplitude");
	g_peaks_time->GetXaxis()->SetTitle("Time");
	g_peaks_time->GetYaxis()->SetTitle("Kick amplitude [V]");
	g_peaks_time->GetXaxis()->SetTimeFormat("%H:%M");
	g_peaks_time->GetXaxis()->SetTimeOffset(0);
	g_peaks_time->GetXaxis()->SetTimeDisplay(1);
	g_peaks_time->GetYaxis()->SetRangeUser(0.6,0.9);
	g_peaks_time->Draw("APLZ");
	gPad->SetGridy();

	can_amplitudes->cd(2);
	for (int i=0; i<Nrepetitions; i++){
		g_peaks[i]->SetTitle("Kick amplitude");
		g_peaks[i]->GetXaxis()->SetTitle("Kick number");
		g_peaks[i]->GetYaxis()->SetTitle("Kick amplitude [V]");
		g_peaks[i]->GetYaxis()->SetRangeUser(0.6,0.9);
		g_peaks[i]->SetMarkerStyle(20);
		g_peaks[i]->SetMarkerColor(i+1);
		g_peaks[i]->Draw(i==0?"APLZ":"PLZ");
	}
	gPad->SetGridy();


    TCanvas* can_amplitudes2 = new TCanvas("can_amplitudes2","",1800,600);
	g_peaks_time_normalized->SetMarkerStyle(20);
	g_peaks_time_normalized->SetTitle("Kick amplitude");
	g_peaks_time_normalized->GetXaxis()->SetTitle("Time");
	g_peaks_time_normalized->GetYaxis()->SetTitle("Kick amplitude [arb. u.]");
	g_peaks_time_normalized->GetXaxis()->SetTimeFormat("%H:%M");
	g_peaks_time_normalized->GetXaxis()->SetTimeOffset(0);
	g_peaks_time_normalized->GetXaxis()->SetTimeDisplay(1);
	g_peaks_time_normalized->GetYaxis()->SetRangeUser(0.7,1.1);
	g_peaks_time_normalized->Draw("APLZ");
	gPad->SetGridy();

    TCanvas* can_kicksAB = new TCanvas("can_kicksAB","",800,600);
    for (int i=0; i<Nfiles; i++){
    	g_kicksA[i]->SetLineColor(kBlack);
    	g_kicksA[i]->SetTitle("All kicks");
    	g_kicksA[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicksA[i]->GetYaxis()->SetTitle("Trace [V]");
    	g_kicksA[i]->GetXaxis()->SetRangeUser(0,6);
    	g_kicksA[i]->GetYaxis()->SetRangeUser(-0.6,0.6);
    	g_kicksA[i]->Draw(i==0?"APL":"PL");
    }
    for (int i=0; i<Nfiles; i++){
    	g_kicksB[i]->SetLineColor(kBlue);
    	g_kicksB[i]->SetTitle("All kicks");
    	g_kicksB[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicksB[i]->GetYaxis()->SetTitle("Trace [V]");
    	g_kicksB[i]->GetXaxis()->SetRangeUser(0,6);
    	g_kicksB[i]->Draw("PL");
    }

	new TCanvas();
	g_A_time->SetTitle("Peak A");
	g_B_time->SetTitle("Peak B");
	g_AB_time->SetTitle("A+B");
	g_A_time->SetLineWidth(2);
	g_B_time->SetLineWidth(2);
	g_AB_time->SetLineWidth(2);
	g_A_time->SetLineWidth(2);
	g_B_time->SetLineWidth(2);
	g_AB_time->SetLineWidth(2);
	g_A_time->SetLineColor(1);
	g_B_time->SetLineColor(2);
	g_AB_time->SetLineColor(3);
	g_A_time->GetXaxis()->SetTimeFormat("%H:%M");
	g_A_time->GetXaxis()->SetTimeOffset(0);
	g_A_time->GetXaxis()->SetTimeDisplay(1);
	g_B_time->GetXaxis()->SetTimeFormat("%H:%M");
	g_B_time->GetXaxis()->SetTimeOffset(0);
	g_B_time->GetXaxis()->SetTimeDisplay(1);
	g_AB_time->GetXaxis()->SetTimeFormat("%H:%M");
	g_AB_time->GetXaxis()->SetTimeOffset(0);
	g_AB_time->GetXaxis()->SetTimeDisplay(1);
	g_A_time->GetYaxis()->SetRangeUser(0,1);
	g_A_time->Draw("APL");
	g_B_time->Draw("PL");
	g_AB_time->Draw("PL");
	gPad->BuildLegend();


	new TCanvas();
	g_kicks[0]->Draw("AL");
	g_kicks_linearized[0]->Draw("L");








	if (saveOutput){
		for (int i=0; i<8; i++){
    		g_kicks_normalized[i]->Write(Form("normalized_kick_%d",i+1));
    	}
		for (int i=0; i<Nfiles; i++){
    		g_kicks_aligned[i]->Write(Form("kick%d_rep%d",1+i%8,1+(i-i%8)/8));
    	}
    	g_peaks_time->Write("kick_amplitude_trend");
		fout->Write();
		fout->Close();
	}
}
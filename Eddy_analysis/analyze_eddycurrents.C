#include "../analysis_tools.C"


void analyze_eddycurrents(TString folder, TString output_file, int Nfilesmax = -1){

	/*
	if (ROOT_VERSION_CODE < ROOT_VERSION(6,24,0)){
		cout<<"ERROR! You are using a ROOT version older than 6.24\n";
		cout<<"If you are on gm2ita, please execute:\n\n";
		cout<<"  source /opt/cernroot/bin/thisroot.sh\n\n";
		return;
	}
	*/
	
	double firstkick_max = 20.0;
	double kick_dt_guess = 10.0;
	double kick_trigger = -200.0;
	double polarity = 1;

	//Begin reading of files
	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	float tstart = 0.;
	float tend = 100.;
	float dt = (tend-tstart)/Nlines;
	float t_before = -1.0;
	float t_after = 8.0;
	int Nkickbins = (t_after-t_before)/dt;
	float t_after8 = 22.0;
	int Nkickbins8 = (t_after8-t_before)/dt;

	double SNR_th1 = 10.0;
	double SNR_th2 = 15.0;

	TFile* fout = new TFile(output_file,"recreate");

	TH2F* h2_fulltrace = new TH2F("h2_fulltrace","Trace;Time [ms];Voltage [mV]",Nlines,tstart,tend,200,-100,100);
	TProfile* g_fulltrace = new TProfile("trace","Trace;Time [ms];Voltage [mV]",Nlines,tstart,tend);
	TProfile** g_trace_kicks = new TProfile*[8];
	for (int i=0; i<8; i++){
		TString hname = Form("trace_kick%d",i+1);
		TString htitle = Form("Trace Kick %d;Time [ms];Voltage [mV]",i+1);
		g_trace_kicks[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
	}
	TProfile* g_trace_kick8long = new TProfile("trace_kick8long","Trace Kick 8;Time [ms];Voltage [mV]",Nkickbins8,t_before,t_after8);

	TH2F* h2_fulltrace_SNR2 = new TH2F("h2_fulltrace_SNR2",Form("Trace (SNR > %.1f);Time [ms];Voltage [mV]",SNR_th2),Nlines,tstart,tend,200,-100,100);
	TProfile* g_fulltrace_SNR0 = new TProfile("trace_SNR0",Form("Trace (SNR < %.1f);Time [ms];Voltage [mV]",SNR_th1),Nlines,tstart,tend);
	TProfile* g_fulltrace_SNR1 = new TProfile("trace_SNR1",Form("Trace (SNR > %.1f);Time [ms];Voltage [mV]",SNR_th1),Nlines,tstart,tend);
	TProfile* g_fulltrace_SNR2 = new TProfile("trace_SNR2",Form("Trace (SNR > %.1f);Time [ms];Voltage [mV]",SNR_th2),Nlines,tstart,tend);
	TProfile** g_trace_kicks_SNR0 = new TProfile*[8];
	TProfile** g_trace_kicks_SNR1 = new TProfile*[8];
	TProfile** g_trace_kicks_SNR2 = new TProfile*[8];
	TProfile* g_fulltrace_SNR2_timealigned = new TProfile("trace_SNR2_timealigned",Form("Aligned trace (SNR > %.1f);Time [ms];Voltage [mV]",SNR_th2),Nlines,tstart-5,tend-5);
	
	TProfile* g_fulltrace_skipped = new TProfile("trace_skipped","Skipped traces;Time [ms];Voltage [mV]",Nlines,tstart-5,tend-5);


	for (int i=0; i<8; i++){
		TString hname = Form("trace_SNR0_kick%d",i+1);
		TString htitle = Form("Trace (SNR < %.1f) Kick %d;Time [ms];Voltage [mV]",SNR_th1,i+1);
		g_trace_kicks_SNR0[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
		
		hname = Form("trace_SNR1_kick%d",i+1);
		htitle = Form("Trace (SNR > %.1f) Kick %d;Time [ms];Voltage [mV]",SNR_th1,i+1);
		g_trace_kicks_SNR1[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
		
		hname = Form("trace_SNR2_kick%d",i+1);
		htitle = Form("Trace (SNR > %.1f) Kick %d;Time [ms];Voltage [mV]",SNR_th2,i+1);
		g_trace_kicks_SNR2[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
	}
	TProfile* g_trace_kick8long_SNR0 = new TProfile("trace_SNR0_kick8long",Form("Trace (SNR < %.1f) Kick 8;Time [ms];Voltage [mV]",SNR_th1),Nkickbins8,t_before,t_after8);
	TProfile* g_trace_kick8long_SNR1 = new TProfile("trace_SNR1_kick8long",Form("Trace (SNR > %.1f) Kick 8;Time [ms];Voltage [mV]",SNR_th1),Nkickbins8,t_before,t_after8);
	TProfile* g_trace_kick8long_SNR2 = new TProfile("trace_SNR2_kick8long",Form("Trace (SNR > %.1f) Kick 8;Time [ms];Voltage [mV]",SNR_th2),Nkickbins8,t_before,t_after8);



	
	TProfile** g_trace_trend = new TProfile*[10];
	for (int i=0; i<10; i++){
		TString hname = Form("trace_trend_%d_%d",i*10,(i+1)*10);
		TString htitle = Form("Trace [%d-%d]%%;Time [ms];Voltage [mV]",i*10,(i+1)*10);
		g_trace_trend[i] = new TProfile(hname,htitle,Nlines,tstart,tend);
	}

	TProfile* g_trace_kickavg = (TProfile*)g_trace_kicks[0]->Clone("trace_avg");
	g_trace_kickavg->Reset();
	g_trace_kickavg->SetTitle("Kick average");

	TGraph* g_trend_A = new TGraph();
	TGraph* g_trend_B = new TGraph();
	TGraph* g_trend_ABdiff = new TGraph();
	TGraphErrors* g_trend_baseline = new TGraphErrors();
	TGraph* g_trend_stddev = new TGraph();
	TGraphErrors* g_trend_blumlein = new TGraphErrors();
	TGraph* g_trend_SNR = new TGraph();

	TGraph* g_correlation_AB_blumlein = new TGraph();
	TGraph* g_correlation_ABdiff_blumlein = new TGraph();
	TGraph* g_correlation_AB_SNR = new TGraph();
	TGraph* g_correlation_ABdiff_SNR = new TGraph();

	TGraph* g_blumlein_temp = new TGraph();
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.5,0.0);
	double blumlein_fit_start = -0.4;
	double blumlein_fit_end = -0.2; 

	for (int fi=0; fi<Nfiles; fi++){

		if (Nfilesmax > 0 && fi >= Nfilesmax) break;

		TString fname = files[fi];
		//if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";
		cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";
		vector<pair<double,double>> trace;

		TDatime datetime = getFileTime(fname);

		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		//Reading file
		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_A = traces[map_varnames["Channel A"]];
		vector<double> trace_B = traces[map_varnames["Channel B"]];
		vector<double> trace_C = traces[map_varnames["Channel C"]];
		vector<double> trace_avgC = traces[map_varnames["average(C)"]];
	
		if (trace_time.size() != Nlines){
			cout<<"Warning! File "<<fname<<" has "<<trace_time.size()<<" lines instead of the expected "<<Nlines<<"!\n";
		}

		//Find kick polarity
		for (int i=0; i<trace_avgC.size(); i++){
			if (trace_time[i] > firstkick_max) break;
			if (trace_avgC[i] < -2000){
				kick_trigger = -abs(kick_trigger);
				polarity = 1;
				break;
			}
			if (trace_avgC[i] > 2000){
				kick_trigger = abs(kick_trigger);
				polarity = -1;
				break;
			}
		}

		//Find first kick
		double first_kick_guess = -1;
		for (int i=1; i<trace_time.size(); i++){
			if ((polarity>0 && trace_avgC[i-1] > kick_trigger && trace_avgC[i] < kick_trigger)||
				(polarity<0 && trace_avgC[i-1] < kick_trigger && trace_avgC[i] > kick_trigger)){
				first_kick_guess = trace_time[i];
				break;
			}
		}

		bool skipped = false;
		if (first_kick_guess < 0){
			cout<<"Warning! Could not find first kick in file "<<fname<<". Skipping...\n";
			skipped = true;
		} else if (first_kick_guess > firstkick_max){
			cout<<"Warning! First kick found at "<<first_kick_guess<<" in file "<<fname<<". Skipping...\n";
			skipped = true;
		}

		if (skipped){
			//subtract baseline anyway from 0 to 5 ms
			double avgC_skipbaseline = 0;
			int Nlines_skipbaseline = 0;
			for (int i=0; i<trace_time.size(); i++){
				if (trace_time[i] < 5){
					avgC_skipbaseline += trace_avgC[i];
					Nlines_skipbaseline++;
				}
			}
			avgC_skipbaseline /= Nlines_skipbaseline;
			
			for (int i=0; i<trace_time.size(); i++){
				g_fulltrace_skipped->Fill(trace_time[i],trace_avgC[i] - avgC_skipbaseline);
			}

			continue;
		}

		double A_avg = 0;
		double B_avg = 0;
		double avgC_baseline = 0;
		double avgC_stddev = 0;
		double avgC_squared = 0;
		double avgC_baselineerr = 0;
		int Nlines_baseline = 0;

		//Measure A, B, baseline, blumlein
		g_blumlein_temp->Set(0);
		for (int i=0; i<trace_time.size(); i++){

			A_avg += trace_A[i];
			B_avg += trace_B[i];

			if (trace_time[i] < first_kick_guess-0.55){
				avgC_baseline += trace_avgC[i];
				avgC_squared += trace_avgC[i]*trace_avgC[i];
				Nlines_baseline++;
			}

			if (trace_time[i] > first_kick_guess-0.55 && trace_time[i] < first_kick_guess){
				g_blumlein_temp->SetPoint(g_blumlein_temp->GetN(),trace_time[i]-first_kick_guess,trace_avgC[i]);
			}
		}

		A_avg /= trace_time.size();
		B_avg /= trace_time.size();
		avgC_baseline /= Nlines_baseline;
		avgC_squared /= Nlines_baseline;
		avgC_stddev = sqrt(avgC_squared - avgC_baseline*avgC_baseline);
		avgC_baselineerr = avgC_stddev / sqrt(Nlines_baseline);

		f_blumlein->SetParameters(60.0,0.5*(blumlein_fit_start+blumlein_fit_end),-1000.0);
		TFitResultPtr fit_blumlein = g_blumlein_temp->Fit("f_blumlein","QS","",blumlein_fit_start,blumlein_fit_end);
		double peak = avgC_baseline;
		double peakerr = 0;
		if (fit_blumlein>=0){
			peak = fit_blumlein->Parameter(0);
			peakerr = fit_blumlein->ParError(0);
		}
		double blumlein = peak - avgC_baseline;
		double blumerr = sqrt(peakerr*peakerr + avgC_baselineerr*avgC_baselineerr);
		double AB = A_avg + B_avg;
		double ABdiff = B_avg - A_avg;
		double SNR = blumlein/avgC_stddev;

		g_trend_A->SetPoint(g_trend_A->GetN(),datetime.Convert(),A_avg);
		g_trend_B->SetPoint(g_trend_B->GetN(),datetime.Convert(),B_avg);
		g_trend_ABdiff->SetPoint(g_trend_ABdiff->GetN(),datetime.Convert(),ABdiff);
		g_trend_baseline->SetPoint(g_trend_baseline->GetN(),datetime.Convert(),avgC_baseline);
		g_trend_baseline->SetPointError(g_trend_baseline->GetN()-1,0,avgC_baselineerr);
		g_trend_stddev->SetPoint(g_trend_stddev->GetN(),datetime.Convert(),avgC_stddev);
		g_trend_blumlein->SetPoint(g_trend_blumlein->GetN(),datetime.Convert(),blumlein);
		g_trend_blumlein->SetPointError(g_trend_blumlein->GetN()-1,0,blumerr);
		g_trend_SNR->SetPoint(g_trend_SNR->GetN(),datetime.Convert(),SNR);

		g_correlation_AB_blumlein->SetPoint(g_correlation_AB_blumlein->GetN(),AB,blumlein);
		g_correlation_ABdiff_blumlein->SetPoint(g_correlation_ABdiff_blumlein->GetN(),ABdiff,blumlein);
		g_correlation_AB_SNR->SetPoint(g_correlation_AB_SNR->GetN(),AB,SNR);
		g_correlation_ABdiff_SNR->SetPoint(g_correlation_ABdiff_SNR->GetN(),ABdiff,SNR);

		///////////////////////////////////////
		//Correct the trace (subtract baseline)
		for (int i=0; i<trace_avgC.size(); i++){
			trace_avgC[i] -= avgC_baseline;
		}

		//Fill the full trace
		for (int i=0; i<trace_time.size(); i++){
			h2_fulltrace->Fill(trace_time[i],trace_avgC[i]);
			g_fulltrace->Fill(trace_time[i],trace_avgC[i]);
			if (SNR < SNR_th1) g_fulltrace_SNR0->Fill(trace_time[i],trace_avgC[i]);
			if (SNR > SNR_th1) g_fulltrace_SNR1->Fill(trace_time[i],trace_avgC[i]);
			if (SNR > SNR_th2){
				h2_fulltrace_SNR2->Fill(trace_time[i],trace_avgC[i]);
				g_fulltrace_SNR2->Fill(trace_time[i],trace_avgC[i]);
			}
			int Nfi_10p = Nfiles / 10;
			int fi_p = (fi-fi%Nfi_10p)/Nfi_10p;
			if (fi_p<10){
				g_trace_trend[fi_p]->Fill(trace_time[i],trace_avgC[i]);
			}
		}

		//Fill individual kicks
		vector<double> kick_timings;
		for (int i=0; i<8; i++){
			double this_guess = first_kick_guess+i*kick_dt_guess-0.4;

			for (int j=0; j<trace_time.size(); j++){
				if (trace_time[j] > this_guess){
					if ((polarity>0 && trace_avgC[j-1] > kick_trigger && trace_avgC[j] < kick_trigger)||
						(polarity<0 && trace_avgC[j-1] < kick_trigger && trace_avgC[j] > kick_trigger)){
						double fraction = (kick_trigger-trace_avgC[j-1])/(trace_avgC[j]-trace_avgC[j-1]);
						double kick_time = trace_time[j-1] + dt*fraction;
						kick_timings.push_back(kick_time);
						for (int k=0; k<trace_time.size(); k++){
							if (trace_time[k] >= kick_time + t_before && trace_time[k] < kick_time + t_after){
								g_trace_kicks[i]->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR < SNR_th1) g_trace_kicks_SNR0[i]->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR > SNR_th1) g_trace_kicks_SNR1[i]->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR > SNR_th2) g_trace_kicks_SNR2[i]->Fill(trace_time[k]-kick_time,trace_avgC[k]);
							}
							if (i==7 && trace_time[k] >= kick_time + t_before && trace_time[k] < kick_time + t_after8){
								g_trace_kick8long->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR < SNR_th1) g_trace_kick8long_SNR0->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR > SNR_th1) g_trace_kick8long_SNR1->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								if (SNR > SNR_th2) g_trace_kick8long_SNR2->Fill(trace_time[k]-kick_time,trace_avgC[k]);
							}
						}
						break;
					}
				}
			}
		}

		//Fill the aligned trace
		for (int i=0; i<trace_time.size(); i++){
			if (SNR > SNR_th2){
				g_fulltrace_SNR2_timealigned->Fill(trace_time[i]-kick_timings[0],trace_avgC[i]);
			}
		}

	}

	for (int i=0; i<8; i++){
		g_trace_kickavg->Add(g_trace_kicks[i]);
	}


	g_trend_A->SetName("g_trend_A");
	g_trend_B->SetName("g_trend_B");
	g_trend_ABdiff->SetName("g_trend_ABdiff");
	g_trend_baseline->SetName("g_trend_baseline");
	g_trend_stddev->SetName("g_trend_stddev");
	g_trend_blumlein->SetName("g_trend_blumlein");
	g_trend_SNR->SetName("g_trend_SNR");
	g_correlation_AB_blumlein->SetName("g_correlation_AB_blumlein");
	g_correlation_ABdiff_blumlein->SetName("g_correlation_ABdiff_blumlein");
	g_correlation_AB_SNR->SetName("g_correlation_AB_SNR");
	g_correlation_ABdiff_SNR->SetName("g_correlation_ABdiff_SNR");

	g_trend_A->SetTitle("Trend of channel A;Time;A [V]");
	g_trend_B->SetTitle("Trend of channel B;Time;B [V]");
	g_trend_ABdiff->SetTitle("Trend of B-A;Time;B-A [V]");
	g_trend_baseline->SetTitle("Trend of avgC baseline;Time;Baseline [mV]");
	g_trend_stddev->SetTitle("Trend of avgC StdDev;Time;StdDev [mV]");
	g_trend_blumlein->SetTitle("Trend of blumlein peak;Time;Blumlein [mV]");
	g_trend_SNR->SetTitle("Trend of SNR;Time;Blumlein / StdDev");
	g_correlation_AB_blumlein->SetTitle("Correlation blumlein vs A+B;A+B [V];Blumlein [mV]");
	g_correlation_ABdiff_blumlein->SetTitle("Correlation blumlein vs B-A;B-A [V];Blumlein [mV]");
	g_correlation_AB_SNR->SetTitle("Correlation SNR vs A+B;A+B [V];Blumlein SNR");
	g_correlation_ABdiff_SNR->SetTitle("Correlation SNR vs B-A;B-A [V];Blumlein SNR");

	g_trend_A->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_A->GetXaxis()->SetTimeOffset(0);
	g_trend_A->GetXaxis()->SetTimeDisplay(1);
	g_trend_B->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_B->GetXaxis()->SetTimeOffset(0);
	g_trend_B->GetXaxis()->SetTimeDisplay(1);
	g_trend_ABdiff->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_ABdiff->GetXaxis()->SetTimeOffset(0);
	g_trend_ABdiff->GetXaxis()->SetTimeDisplay(1);
	g_trend_baseline->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_baseline->GetXaxis()->SetTimeOffset(0);
	g_trend_baseline->GetXaxis()->SetTimeDisplay(1);
	g_trend_stddev->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_stddev->GetXaxis()->SetTimeOffset(0);
	g_trend_stddev->GetXaxis()->SetTimeDisplay(1);
	g_trend_blumlein->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_blumlein->GetXaxis()->SetTimeOffset(0);
	g_trend_blumlein->GetXaxis()->SetTimeDisplay(1);
	g_trend_SNR->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_SNR->GetXaxis()->SetTimeOffset(0);
	g_trend_SNR->GetXaxis()->SetTimeDisplay(1);

	g_trend_A->Write();
	g_trend_B->Write();
	g_trend_ABdiff->Write();
	g_trend_baseline->Write();
	g_trend_stddev->Write();
	g_trend_blumlein->Write();
	g_trend_SNR->Write();

	g_correlation_AB_blumlein->Write();
	g_correlation_ABdiff_blumlein->Write();
	g_correlation_AB_SNR->Write();
	g_correlation_ABdiff_SNR->Write();

	fout->Write();
	fout->Close();
	cout<<"\nDone! Created output file "<<output_file<<"\n";
	return;
}

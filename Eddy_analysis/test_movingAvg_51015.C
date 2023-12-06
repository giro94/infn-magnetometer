#include "../analysis_tools.C"



void test_movingAvg_51015(TString input_file, TString output_file=""){

	TFile* f = TFile::Open(input_file);


	bool saveOutput = false;
	if (output_file != ""){
		saveOutput = true;
	}

	gStyle->SetOptStat(0);

	TH1D* trace = ((TProfile*)f->Get("trace_SNR2"))->ProjectionX();
	TH1D** kicks = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks[i] = ((TProfile*)f->Get(Form("trace_SNR2_kick%d",i+1)))->ProjectionX();
	}
	TH1D* kick8long = ((TProfile*)f->Get("trace_SNR2_kick8long"))->ProjectionX();

	for (int bn=1; bn<=trace->GetNbinsX(); bn++){
		if (trace->GetBinContent(bn) < -50){
			trace->SetBinContent(bn,0);
		}
	}

	for (int i=0; i<8; i++){
		for (int bn=1; bn<=kicks[i]->GetNbinsX(); bn++){
			if (kicks[i]->GetBinContent(bn) < -50){
				kicks[i]->SetBinContent(bn,0);
			}
		}
	}

	for (int bn=1; bn<=kick8long->GetNbinsX(); bn++){
		if (kick8long->GetBinContent(bn) < -50){
			kick8long->SetBinContent(bn,0);
		}
	}

	int running_width;

	running_width = 5;
	TH1D* trace_ra5 = runningAverage(trace,running_width,true);
	trace_ra5->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),running_width));	
	TH1D** kicks_ra5 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra5[i] = runningAverage(kicks[i],running_width,true);
		kicks_ra5[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),running_width));	
	}
	TH1D* kick8long_ra5 = runningAverage(kick8long,running_width,true);
	kick8long_ra5->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),running_width));	

	running_width = 10;
	TH1D* trace_ra10 = runningAverage(trace,running_width,true);
	trace_ra10->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),running_width));	
	TH1D** kicks_ra10 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra10[i] = runningAverage(kicks[i],running_width,true);
		kicks_ra10[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),running_width));	
	}
	TH1D* kick8long_ra10 = runningAverage(kick8long,running_width,true);
	kick8long_ra10->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),running_width));	


	running_width = 15;
	TH1D* trace_ra15 = runningAverage(trace,running_width,true);
	trace_ra15->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),running_width));	
	TH1D** kicks_ra15 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra15[i] = runningAverage(kicks[i],running_width,true);
		kicks_ra15[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),running_width));	
	}
	TH1D* kick8long_ra15 = runningAverage(kick8long,running_width,true);
	kick8long_ra15->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),running_width));	

	TH1D* trace_ra_5_10_15 = runningAverage(runningAverage(runningAverage(trace,5,true),10,true),15,true);
	trace_ra_5_10_15->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),51015));	
	TH1D** kicks_ra_5_10_15 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra_5_10_15[i] = runningAverage(runningAverage(runningAverage(kicks[i],5,true),10,true),15,true);
		kicks_ra_5_10_15[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),51015));	
	}
	TH1D* kick8long_ra_5_10_15 = runningAverage(runningAverage(runningAverage(kick8long,5,true),10,true),15,true);
	kick8long_ra_5_10_15->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),51015));	

	new TCanvas();
	trace->SetLineWidth(2);
	trace_ra5->SetLineWidth(2);
	trace_ra10->SetLineWidth(2);
	trace_ra15->SetLineWidth(2);
	trace_ra_5_10_15->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_ra5->SetLineColor(kRed);
	trace_ra10->SetLineColor(kOrange);
	trace_ra15->SetLineColor(kGreen+2);
	trace_ra_5_10_15->SetLineColor(kViolet);
	trace->Draw("HIST");
	trace_ra5->Draw("HIST SAME");
	trace_ra10->Draw("HIST SAME");
	trace_ra15->Draw("HIST SAME");
	trace_ra_5_10_15->Draw("HIST SAME");


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra5[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra10[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra15[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra_5_10_15[i]->GetYaxis()->SetRangeUser(-60,60);
	}

	TCanvas* can_kicks = new TCanvas("can_kicks","",1600,800);
	can_kicks->Divide(2,1);
	can_kicks->cd(1);
	kicks[0]->SetLineWidth(2);
	kicks_ra5[0]->SetLineWidth(2);
	kicks_ra10[0]->SetLineWidth(2);
	kicks_ra15[0]->SetLineWidth(2);
	kicks_ra_5_10_15[0]->SetLineWidth(2);
	kicks[0]->SetLineColor(kBlue);
	kicks_ra5[0]->SetLineColor(kRed);
	kicks_ra10[0]->SetLineColor(kOrange);
	kicks_ra15[0]->SetLineColor(kGreen+2);
	kicks_ra_5_10_15[0]->SetLineColor(kViolet);
	kicks[0]->Draw("HIST");
	kicks_ra5[0]->Draw("HIST SAME");
	kicks_ra10[0]->Draw("HIST SAME");
	kicks_ra15[0]->Draw("HIST SAME");
	kicks_ra_5_10_15[0]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks->cd(2);
	kicks[7]->SetLineWidth(2);
	kicks_ra5[7]->SetLineWidth(2);
	kicks_ra10[7]->SetLineWidth(2);
	kicks_ra15[7]->SetLineWidth(2);
	kicks_ra_5_10_15[7]->SetLineWidth(2);
	kicks[7]->SetLineColor(kBlue);
	kicks_ra5[7]->SetLineColor(kRed);
	kicks_ra10[7]->SetLineColor(kOrange);
	kicks_ra15[7]->SetLineColor(kGreen+2);
	kicks_ra_5_10_15[7]->SetLineColor(kViolet);
	kicks[7]->Draw("HIST");
	kicks_ra5[7]->Draw("HIST SAME");
	kicks_ra10[7]->Draw("HIST SAME");
	kicks_ra15[7]->Draw("HIST SAME");
	kicks_ra_5_10_15[7]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);


	TH1D* kick1_zoom = (TH1D*)kicks[0]->Clone("kick1_zoom");
	TH1D* kick1_zoom_ra5 = (TH1D*)kicks_ra5[0]->Clone("kick1_zoom_ra5");
	TH1D* kick1_zoom_ra10 = (TH1D*)kicks_ra10[0]->Clone("kick1_zoom_ra10");
	TH1D* kick1_zoom_ra15 = (TH1D*)kicks_ra15[0]->Clone("kick1_zoom_ra15");
	TH1D* kick1_zoom_ra_5_10_15 = (TH1D*)kicks_ra_5_10_15[0]->Clone("kick1_zoom_ra_5_10_15");

	TH1D* kick8_zoom = (TH1D*)kicks[7]->Clone("kick8_zoom");
	TH1D* kick8_zoom_ra5 = (TH1D*)kicks_ra5[7]->Clone("kick8_zoom_ra5");
	TH1D* kick8_zoom_ra10 = (TH1D*)kicks_ra10[7]->Clone("kick8_zoom_ra10");
	TH1D* kick8_zoom_ra15 = (TH1D*)kicks_ra15[7]->Clone("kick8_zoom_ra15");
	TH1D* kick8_zoom_ra_5_10_15 = (TH1D*)kicks_ra_5_10_15[7]->Clone("kick8_zoom_ra_5_10_15");

	double trace_offset = 25;
	for (int bn=1; bn<=kick1_zoom->GetNbinsX(); bn++){
		kick1_zoom_ra5->SetBinContent(bn,kick1_zoom_ra5->GetBinContent(bn) - trace_offset);
		kick1_zoom_ra10->SetBinContent(bn,kick1_zoom_ra10->GetBinContent(bn) - 2*trace_offset);
		kick1_zoom_ra15->SetBinContent(bn,kick1_zoom_ra15->GetBinContent(bn) - 3*trace_offset);
		kick1_zoom_ra_5_10_15->SetBinContent(bn,kick1_zoom_ra_5_10_15->GetBinContent(bn) - 4*trace_offset);

		kick8_zoom_ra5->SetBinContent(bn,kick8_zoom_ra5->GetBinContent(bn) - trace_offset);
		kick8_zoom_ra10->SetBinContent(bn,kick8_zoom_ra10->GetBinContent(bn) - 2*trace_offset);
		kick8_zoom_ra15->SetBinContent(bn,kick8_zoom_ra15->GetBinContent(bn) - 3*trace_offset);
		kick8_zoom_ra_5_10_15->SetBinContent(bn,kick8_zoom_ra_5_10_15->GetBinContent(bn) - 4*trace_offset);
	}

	double x_zoom_from = 0.0;
	double x_zoom_to = 1.0;
	kick1_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra5->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra10->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra15->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra_5_10_15->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra5->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra10->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra15->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra_5_10_15->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);

	kick1_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra5->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra10->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra15->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra_5_10_15->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra5->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra10->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra15->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra_5_10_15->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);

	TCanvas* can_kicks_zoom = new TCanvas("can_kicks_zoom","",1800,900);
	can_kicks_zoom->Divide(2,1);
	can_kicks_zoom->cd(1);
	kick1_zoom->Draw("HIST");
	kick1_zoom_ra5->Draw("HIST SAME");
	kick1_zoom_ra10->Draw("HIST SAME");
	kick1_zoom_ra15->Draw("HIST SAME");
	kick1_zoom_ra_5_10_15->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks_zoom->cd(2);
	kick8_zoom->Draw("HIST");
	kick8_zoom_ra5->Draw("HIST SAME");
	kick8_zoom_ra10->Draw("HIST SAME");
	kick8_zoom_ra15->Draw("HIST SAME");
	kick8_zoom_ra_5_10_15->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra5[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra10[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra15[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra_5_10_15[i]->GetYaxis()->SetRangeUser(-60,60);
	}


	TGraphErrors* g_off = new TGraphErrors();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	TGraphErrors* g_off_ra5 = new TGraphErrors();
	TGraphErrors* g_amp_ra5 = new TGraphErrors();
	TGraphErrors* g_tau_ra5 = new TGraphErrors();
	TGraphErrors* g_off_ra10 = new TGraphErrors();
	TGraphErrors* g_amp_ra10 = new TGraphErrors();
	TGraphErrors* g_tau_ra10 = new TGraphErrors();
	TGraphErrors* g_off_ra15 = new TGraphErrors();
	TGraphErrors* g_amp_ra15 = new TGraphErrors();
	TGraphErrors* g_tau_ra15 = new TGraphErrors();
	TGraphErrors* g_off_ra_5_10_15 = new TGraphErrors();
	TGraphErrors* g_amp_ra_5_10_15 = new TGraphErrors();
	TGraphErrors* g_tau_ra_5_10_15 = new TGraphErrors();
	g_off->SetTitle("Offset;Kick;Offset [mV]");
	g_amp->SetTitle("Amplitude;Kick;Amplitude [mV]");
	g_tau->SetTitle("Lifetime;Kick;Lifetime [#mus]");
	g_off_ra5->SetTitle("Offset (RA 10);Kick;Offset [mV]");
	g_amp_ra5->SetTitle("Amplitude (RA 10);Kick;Amplitude [mV]");
	g_tau_ra5->SetTitle("Lifetime (RA 10);Kick;Lifetime [#mus]");
	g_off_ra10->SetTitle("Offset (RA 20);Kick;Offset [mV]");
	g_amp_ra10->SetTitle("Amplitude (RA 20);Kick;Amplitude [mV]");
	g_tau_ra10->SetTitle("Lifetime (RA 20);Kick;Lifetime [#mus]");
	g_off_ra15->SetTitle("Offset (RA 30);Kick;Offset [mV]");
	g_amp_ra15->SetTitle("Amplitude (RA 30);Kick;Amplitude [mV]");
	g_tau_ra15->SetTitle("Lifetime (RA 30);Kick;Lifetime [#mus]");
	g_off_ra_5_10_15->SetTitle("Offset (RA 51015);Kick;Offset [mV]");
	g_amp_ra_5_10_15->SetTitle("Amplitude (RA 51015);Kick;Amplitude [mV]");
	g_tau_ra_5_10_15->SetTitle("Lifetime (RA 51015);Kick;Lifetime [#mus]");
	g_off->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_off_ra5->SetMarkerStyle(20);
	g_amp_ra5->SetMarkerStyle(20);
	g_tau_ra5->SetMarkerStyle(20);
	g_off_ra10->SetMarkerStyle(20);
	g_amp_ra10->SetMarkerStyle(20);
	g_tau_ra10->SetMarkerStyle(20);
	g_off_ra15->SetMarkerStyle(20);
	g_amp_ra15->SetMarkerStyle(20);
	g_tau_ra15->SetMarkerStyle(20);
	g_off_ra_5_10_15->SetMarkerStyle(20);
	g_amp_ra_5_10_15->SetMarkerStyle(20);
	g_tau_ra_5_10_15->SetMarkerStyle(20);
	g_off->SetMarkerColor(kBlue);
	g_amp->SetMarkerColor(kBlue);
	g_tau->SetMarkerColor(kBlue);
	g_off_ra5->SetMarkerColor(kRed);
	g_amp_ra5->SetMarkerColor(kRed);
	g_tau_ra5->SetMarkerColor(kRed);
	g_off_ra10->SetMarkerColor(kOrange);
	g_amp_ra10->SetMarkerColor(kOrange);
	g_tau_ra10->SetMarkerColor(kOrange);
	g_off_ra15->SetMarkerColor(kGreen+2);
	g_amp_ra15->SetMarkerColor(kGreen+2);
	g_tau_ra15->SetMarkerColor(kGreen+2);
	g_off_ra_5_10_15->SetMarkerColor(kViolet);
	g_amp_ra_5_10_15->SetMarkerColor(kViolet);
	g_tau_ra_5_10_15->SetMarkerColor(kViolet);


	TF1* f_exp = new TF1("f_exp","[2]-[0]*exp(-x/[1])",0.02,0.7);
	for (int i=0; i<8; i++){
		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick = kicks[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra5 = kicks_ra5[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra5->SetPoint(g_off_ra5->GetN(),i+1,fit_kick_ra5->Parameter(2));
		g_off_ra5->SetPointError(g_off_ra5->GetN()-1,0,fit_kick_ra5->ParError(2));
		g_amp_ra5->SetPoint(g_amp_ra5->GetN(),i+1,fit_kick_ra5->Parameter(0));
		g_amp_ra5->SetPointError(g_amp_ra5->GetN()-1,0,fit_kick_ra5->ParError(0));
		g_tau_ra5->SetPoint(g_tau_ra5->GetN(),i+1,1000*fit_kick_ra5->Parameter(1));
		g_tau_ra5->SetPointError(g_tau_ra5->GetN()-1,0,1000*fit_kick_ra5->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra10 = kicks_ra10[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra10->SetPoint(g_off_ra10->GetN(),i+1,fit_kick_ra10->Parameter(2));
		g_off_ra10->SetPointError(g_off_ra10->GetN()-1,0,fit_kick_ra10->ParError(2));
		g_amp_ra10->SetPoint(g_amp_ra10->GetN(),i+1,fit_kick_ra10->Parameter(0));
		g_amp_ra10->SetPointError(g_amp_ra10->GetN()-1,0,fit_kick_ra10->ParError(0));
		g_tau_ra10->SetPoint(g_tau_ra10->GetN(),i+1,1000*fit_kick_ra10->Parameter(1));
		g_tau_ra10->SetPointError(g_tau_ra10->GetN()-1,0,1000*fit_kick_ra10->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra15 = kicks_ra15[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra15->SetPoint(g_off_ra15->GetN(),i+1,fit_kick_ra15->Parameter(2));
		g_off_ra15->SetPointError(g_off_ra15->GetN()-1,0,fit_kick_ra15->ParError(2));
		g_amp_ra15->SetPoint(g_amp_ra15->GetN(),i+1,fit_kick_ra15->Parameter(0));
		g_amp_ra15->SetPointError(g_amp_ra15->GetN()-1,0,fit_kick_ra15->ParError(0));
		g_tau_ra15->SetPoint(g_tau_ra15->GetN(),i+1,1000*fit_kick_ra15->Parameter(1));
		g_tau_ra15->SetPointError(g_tau_ra15->GetN()-1,0,1000*fit_kick_ra15->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra_5_10_15 = kicks_ra_5_10_15[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra_5_10_15->SetPoint(g_off_ra_5_10_15->GetN(),i+1,fit_kick_ra_5_10_15->Parameter(2));
		g_off_ra_5_10_15->SetPointError(g_off_ra_5_10_15->GetN()-1,0,fit_kick_ra_5_10_15->ParError(2));
		g_amp_ra_5_10_15->SetPoint(g_amp_ra_5_10_15->GetN(),i+1,fit_kick_ra_5_10_15->Parameter(0));
		g_amp_ra_5_10_15->SetPointError(g_amp_ra_5_10_15->GetN()-1,0,fit_kick_ra_5_10_15->ParError(0));
		g_tau_ra_5_10_15->SetPoint(g_tau_ra_5_10_15->GetN(),i+1,1000*fit_kick_ra_5_10_15->Parameter(1));
		g_tau_ra_5_10_15->SetPointError(g_tau_ra_5_10_15->GetN()-1,0,1000*fit_kick_ra_5_10_15->ParError(1));

		//kicks[i]->Add(f_exp,-1);
		//kicks_ra5[i]->Add(f_exp,-1);
		//kicks_ra10[i]->Add(f_exp,-1);
		//kicks_ra15[i]->Add(f_exp,-1);
		//kicks_ra_5_10_15[i]->Add(f_exp,-1);
	}

	TH1D** h1_fft_kicks = new TH1D* [8];
	TH1D** h1_fft_kicks_ra5 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra10 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra15 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra_5_10_15 = new TH1D* [8];
	double fft_xmin = 0.1; //ms from kick
	double fft_xmax = 8.0; 
	double dt = kicks[0]->GetBinWidth(1);
	for (int i=0; i<8; i++){
		TH1 *fft_histogram = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit = SetupFFT(kicks[i], fft_xmin, fft_xmax);
		fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
		h1_fft_kicks[i] = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks[i]->SetStats(0);
		h1_fft_kicks[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks[i]->Scale(1.0 / h1_fft_kicks[i]->Integral());
		h1_fft_kicks[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra5 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra5 = SetupFFT(kicks_ra5[i], fft_xmin, fft_xmax);
		fft_histogram_ra5 = fftResidualInit_ra5->FFT(fft_histogram_ra5,"MAG");
		h1_fft_kicks_ra5[i] = RescaleAxis(fft_histogram_ra5, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra5[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra5[i]->SetStats(0);
		h1_fft_kicks_ra5[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra5[i]->Scale(1.0 / h1_fft_kicks_ra5[i]->Integral());
		h1_fft_kicks_ra5[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra5[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra10 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra10 = SetupFFT(kicks_ra10[i], fft_xmin, fft_xmax);
		fft_histogram_ra10 = fftResidualInit_ra10->FFT(fft_histogram_ra10,"MAG");
		h1_fft_kicks_ra10[i] = RescaleAxis(fft_histogram_ra10, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra10[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra10[i]->SetStats(0);
		h1_fft_kicks_ra10[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra10[i]->Scale(1.0 / h1_fft_kicks_ra10[i]->Integral());
		h1_fft_kicks_ra10[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra10[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra15 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra15 = SetupFFT(kicks_ra15[i], fft_xmin, fft_xmax);
		fft_histogram_ra15 = fftResidualInit_ra15->FFT(fft_histogram_ra15,"MAG");
		h1_fft_kicks_ra15[i] = RescaleAxis(fft_histogram_ra15, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra15[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra15[i]->SetStats(0);
		h1_fft_kicks_ra15[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra15[i]->Scale(1.0 / h1_fft_kicks_ra15[i]->Integral());
		h1_fft_kicks_ra15[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra15[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra_5_10_15 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra_5_10_15 = SetupFFT(kicks_ra_5_10_15[i], fft_xmin, fft_xmax);
		fft_histogram_ra_5_10_15 = fftResidualInit_ra_5_10_15->FFT(fft_histogram_ra_5_10_15,"MAG");
		h1_fft_kicks_ra_5_10_15[i] = RescaleAxis(fft_histogram_ra_5_10_15, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra_5_10_15[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra_5_10_15[i]->SetStats(0);
		h1_fft_kicks_ra_5_10_15[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra_5_10_15[i]->Scale(1.0 / h1_fft_kicks_ra_5_10_15[i]->Integral());
		h1_fft_kicks_ra_5_10_15[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra_5_10_15[i]->GetYaxis()->SetRangeUser(0,0.01);
	}

	TCanvas* can = new TCanvas("can","",1400,800);
	can->Divide(5,2);
	can->cd(1);
	kicks[0]->Draw("HIST");
	can->cd(2);
	kicks_ra5[0]->Draw("HIST");
	can->cd(3);
	kicks_ra10[0]->Draw("HIST");
	can->cd(4);
	kicks_ra15[0]->Draw("HIST");
	can->cd(5);
	kicks_ra_5_10_15[0]->Draw("HIST");
	can->cd(6);
	h1_fft_kicks[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks: "<<h1_fft_kicks[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(7);
	h1_fft_kicks_ra5[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra5: "<<h1_fft_kicks_ra5[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(8);
	h1_fft_kicks_ra10[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra10: "<<h1_fft_kicks_ra10[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(9);
	h1_fft_kicks_ra15[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra15: "<<h1_fft_kicks_ra15[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(10);
	h1_fft_kicks_ra_5_10_15[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra_5_10_15: "<<h1_fft_kicks_ra_5_10_15[0]->Integral()<<"\n";
	gPad->SetLogx();


	TCanvas* can8 = new TCanvas("can8","",1400,800);
	can8->Divide(5,2);
	can8->cd(1);
	kicks[7]->Draw("HIST");
	can8->cd(2);
	kicks_ra5[7]->Draw("HIST");
	can8->cd(3);
	kicks_ra10[7]->Draw("HIST");
	can8->cd(4);
	kicks_ra15[7]->Draw("HIST");
	can8->cd(5);
	kicks_ra_5_10_15[7]->Draw("HIST");
	can8->cd(6);
	h1_fft_kicks[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(7);
	h1_fft_kicks_ra5[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(8);
	h1_fft_kicks_ra10[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(9);
	h1_fft_kicks_ra15[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(10);
	h1_fft_kicks_ra_5_10_15[7]->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can_fits = new TCanvas("can_fits","",1500,500);
	can_fits->Divide(3,1);
	can_fits->cd(1);
	g_off->GetYaxis()->SetRangeUser(-1,1);
	g_off->Draw("APZ");
	g_off_ra5->Draw("PZ");
	g_off_ra10->Draw("PZ");
	g_off_ra15->Draw("PZ");
	g_off_ra_5_10_15->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(2);
	g_amp->GetYaxis()->SetRangeUser(0,25);
	g_amp->Draw("APZ");
	g_amp_ra5->Draw("PZ");
	g_amp_ra10->Draw("PZ");
	g_amp_ra15->Draw("PZ");
	g_amp_ra_5_10_15->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(3);
	g_tau->GetYaxis()->SetRangeUser(20,100);
	g_tau->Draw("APZ");
	g_tau_ra5->Draw("PZ");
	g_tau_ra10->Draw("PZ");
	g_tau_ra15->Draw("PZ");
	g_tau_ra_5_10_15->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);

}

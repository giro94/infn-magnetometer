#include "../analysis_tools.C"


void test_movingAvg_102030(TString input_file, TString output_file=""){

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

	running_width = 20;
	TH1D* trace_ra20 = runningAverage(trace,running_width,true);
	trace_ra20->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),running_width));	
	TH1D** kicks_ra20 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra20[i] = runningAverage(kicks[i],running_width,true);
		kicks_ra20[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),running_width));	
	}
	TH1D* kick8long_ra20 = runningAverage(kick8long,running_width,true);
	kick8long_ra20->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),running_width));	


	running_width = 30;
	TH1D* trace_ra30 = runningAverage(trace,running_width,true);
	trace_ra30->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),running_width));	
	TH1D** kicks_ra30 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra30[i] = runningAverage(kicks[i],running_width,true);
		kicks_ra30[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),running_width));	
	}
	TH1D* kick8long_ra30 = runningAverage(kick8long,running_width,true);
	kick8long_ra30->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),running_width));	

	TH1D* trace_ra_10_20_30 = runningAverage(runningAverage(runningAverage(trace,10,true),20,true),30,true);
	trace_ra_10_20_30->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),102030));	
	TH1D** kicks_ra_10_20_30 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra_10_20_30[i] = runningAverage(runningAverage(runningAverage(kicks[i],10,true),20,true),30,true);
		kicks_ra_10_20_30[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),102030));	
	}
	TH1D* kick8long_ra_10_20_30 = runningAverage(runningAverage(runningAverage(kick8long,10,true),20,true),30,true);
	kick8long_ra_10_20_30->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),102030));	

	new TCanvas();
	trace->SetLineWidth(2);
	trace_ra10->SetLineWidth(2);
	trace_ra20->SetLineWidth(2);
	trace_ra30->SetLineWidth(2);
	trace_ra_10_20_30->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_ra10->SetLineColor(kRed);
	trace_ra20->SetLineColor(kOrange);
	trace_ra30->SetLineColor(kGreen+2);
	trace_ra_10_20_30->SetLineColor(kViolet);
	trace->Draw("HIST");
	trace_ra10->Draw("HIST SAME");
	trace_ra20->Draw("HIST SAME");
	trace_ra30->Draw("HIST SAME");
	trace_ra_10_20_30->Draw("HIST SAME");


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra10[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra20[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra30[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra_10_20_30[i]->GetYaxis()->SetRangeUser(-60,60);
	}

	TCanvas* can_kicks = new TCanvas("can_kicks","",1600,800);
	can_kicks->Divide(2,1);
	can_kicks->cd(1);
	kicks[0]->SetLineWidth(2);
	kicks_ra10[0]->SetLineWidth(2);
	kicks_ra20[0]->SetLineWidth(2);
	kicks_ra30[0]->SetLineWidth(2);
	kicks_ra_10_20_30[0]->SetLineWidth(2);
	kicks[0]->SetLineColor(kBlue);
	kicks_ra10[0]->SetLineColor(kRed);
	kicks_ra20[0]->SetLineColor(kOrange);
	kicks_ra30[0]->SetLineColor(kGreen+2);
	kicks_ra_10_20_30[0]->SetLineColor(kViolet);
	kicks[0]->Draw("HIST");
	kicks_ra10[0]->Draw("HIST SAME");
	kicks_ra20[0]->Draw("HIST SAME");
	kicks_ra30[0]->Draw("HIST SAME");
	kicks_ra_10_20_30[0]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks->cd(2);
	kicks[7]->SetLineWidth(2);
	kicks_ra10[7]->SetLineWidth(2);
	kicks_ra20[7]->SetLineWidth(2);
	kicks_ra30[7]->SetLineWidth(2);
	kicks_ra_10_20_30[7]->SetLineWidth(2);
	kicks[7]->SetLineColor(kBlue);
	kicks_ra10[7]->SetLineColor(kRed);
	kicks_ra20[7]->SetLineColor(kOrange);
	kicks_ra30[7]->SetLineColor(kGreen+2);
	kicks_ra_10_20_30[7]->SetLineColor(kViolet);
	kicks[7]->Draw("HIST");
	kicks_ra10[7]->Draw("HIST SAME");
	kicks_ra20[7]->Draw("HIST SAME");
	kicks_ra30[7]->Draw("HIST SAME");
	kicks_ra_10_20_30[7]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);


	TH1D* kick1_zoom = (TH1D*)kicks[0]->Clone("kick1_zoom");
	TH1D* kick1_zoom_ra10 = (TH1D*)kicks_ra10[0]->Clone("kick1_zoom_ra10");
	TH1D* kick1_zoom_ra20 = (TH1D*)kicks_ra20[0]->Clone("kick1_zoom_ra20");
	TH1D* kick1_zoom_ra30 = (TH1D*)kicks_ra30[0]->Clone("kick1_zoom_ra30");
	TH1D* kick1_zoom_ra_10_20_30 = (TH1D*)kicks_ra_10_20_30[0]->Clone("kick1_zoom_ra_10_20_30");

	TH1D* kick8_zoom = (TH1D*)kicks[7]->Clone("kick8_zoom");
	TH1D* kick8_zoom_ra10 = (TH1D*)kicks_ra10[7]->Clone("kick8_zoom_ra10");
	TH1D* kick8_zoom_ra20 = (TH1D*)kicks_ra20[7]->Clone("kick8_zoom_ra20");
	TH1D* kick8_zoom_ra30 = (TH1D*)kicks_ra30[7]->Clone("kick8_zoom_ra30");
	TH1D* kick8_zoom_ra_10_20_30 = (TH1D*)kicks_ra_10_20_30[7]->Clone("kick8_zoom_ra_10_20_30");

	double trace_offset = 25;
	for (int bn=1; bn<=kick1_zoom->GetNbinsX(); bn++){
		kick1_zoom_ra10->SetBinContent(bn,kick1_zoom_ra10->GetBinContent(bn) - trace_offset);
		kick1_zoom_ra20->SetBinContent(bn,kick1_zoom_ra20->GetBinContent(bn) - 2*trace_offset);
		kick1_zoom_ra30->SetBinContent(bn,kick1_zoom_ra30->GetBinContent(bn) - 3*trace_offset);
		kick1_zoom_ra_10_20_30->SetBinContent(bn,kick1_zoom_ra_10_20_30->GetBinContent(bn) - 4*trace_offset);

		kick8_zoom_ra10->SetBinContent(bn,kick8_zoom_ra10->GetBinContent(bn) - trace_offset);
		kick8_zoom_ra20->SetBinContent(bn,kick8_zoom_ra20->GetBinContent(bn) - 2*trace_offset);
		kick8_zoom_ra30->SetBinContent(bn,kick8_zoom_ra30->GetBinContent(bn) - 3*trace_offset);
		kick8_zoom_ra_10_20_30->SetBinContent(bn,kick8_zoom_ra_10_20_30->GetBinContent(bn) - 4*trace_offset);
	}

	double x_zoom_from = 0.0;
	double x_zoom_to = 1.0;
	kick1_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra10->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra20->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra30->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra_10_20_30->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra10->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra20->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra30->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra_10_20_30->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);

	kick1_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra10->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra20->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra30->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra_10_20_30->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra10->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra20->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra30->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra_10_20_30->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);

	TCanvas* can_kicks_zoom = new TCanvas("can_kicks_zoom","",1800,900);
	can_kicks_zoom->Divide(2,1);
	can_kicks_zoom->cd(1);
	kick1_zoom->Draw("HIST");
	kick1_zoom_ra10->Draw("HIST SAME");
	kick1_zoom_ra20->Draw("HIST SAME");
	kick1_zoom_ra30->Draw("HIST SAME");
	kick1_zoom_ra_10_20_30->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks_zoom->cd(2);
	kick8_zoom->Draw("HIST");
	kick8_zoom_ra10->Draw("HIST SAME");
	kick8_zoom_ra20->Draw("HIST SAME");
	kick8_zoom_ra30->Draw("HIST SAME");
	kick8_zoom_ra_10_20_30->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra10[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra20[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra30[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra_10_20_30[i]->GetYaxis()->SetRangeUser(-60,60);
	}


	TGraphErrors* g_off = new TGraphErrors();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	TGraphErrors* g_off_ra10 = new TGraphErrors();
	TGraphErrors* g_amp_ra10 = new TGraphErrors();
	TGraphErrors* g_tau_ra10 = new TGraphErrors();
	TGraphErrors* g_off_ra20 = new TGraphErrors();
	TGraphErrors* g_amp_ra20 = new TGraphErrors();
	TGraphErrors* g_tau_ra20 = new TGraphErrors();
	TGraphErrors* g_off_ra30 = new TGraphErrors();
	TGraphErrors* g_amp_ra30 = new TGraphErrors();
	TGraphErrors* g_tau_ra30 = new TGraphErrors();
	TGraphErrors* g_off_ra_10_20_30 = new TGraphErrors();
	TGraphErrors* g_amp_ra_10_20_30 = new TGraphErrors();
	TGraphErrors* g_tau_ra_10_20_30 = new TGraphErrors();
	g_off->SetTitle("Offset;Kick;Offset [mV]");
	g_amp->SetTitle("Amplitude;Kick;Amplitude [mV]");
	g_tau->SetTitle("Lifetime;Kick;Lifetime [#mus]");
	g_off_ra10->SetTitle("Offset (RA 10);Kick;Offset [mV]");
	g_amp_ra10->SetTitle("Amplitude (RA 10);Kick;Amplitude [mV]");
	g_tau_ra10->SetTitle("Lifetime (RA 10);Kick;Lifetime [#mus]");
	g_off_ra20->SetTitle("Offset (RA 20);Kick;Offset [mV]");
	g_amp_ra20->SetTitle("Amplitude (RA 20);Kick;Amplitude [mV]");
	g_tau_ra20->SetTitle("Lifetime (RA 20);Kick;Lifetime [#mus]");
	g_off_ra30->SetTitle("Offset (RA 30);Kick;Offset [mV]");
	g_amp_ra30->SetTitle("Amplitude (RA 30);Kick;Amplitude [mV]");
	g_tau_ra30->SetTitle("Lifetime (RA 30);Kick;Lifetime [#mus]");
	g_off_ra_10_20_30->SetTitle("Offset (RA 102030);Kick;Offset [mV]");
	g_amp_ra_10_20_30->SetTitle("Amplitude (RA 102030);Kick;Amplitude [mV]");
	g_tau_ra_10_20_30->SetTitle("Lifetime (RA 102030);Kick;Lifetime [#mus]");
	g_off->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_off_ra10->SetMarkerStyle(20);
	g_amp_ra10->SetMarkerStyle(20);
	g_tau_ra10->SetMarkerStyle(20);
	g_off_ra20->SetMarkerStyle(20);
	g_amp_ra20->SetMarkerStyle(20);
	g_tau_ra20->SetMarkerStyle(20);
	g_off_ra30->SetMarkerStyle(20);
	g_amp_ra30->SetMarkerStyle(20);
	g_tau_ra30->SetMarkerStyle(20);
	g_off_ra_10_20_30->SetMarkerStyle(20);
	g_amp_ra_10_20_30->SetMarkerStyle(20);
	g_tau_ra_10_20_30->SetMarkerStyle(20);
	g_off->SetMarkerColor(kBlue);
	g_amp->SetMarkerColor(kBlue);
	g_tau->SetMarkerColor(kBlue);
	g_off_ra10->SetMarkerColor(kRed);
	g_amp_ra10->SetMarkerColor(kRed);
	g_tau_ra10->SetMarkerColor(kRed);
	g_off_ra20->SetMarkerColor(kOrange);
	g_amp_ra20->SetMarkerColor(kOrange);
	g_tau_ra20->SetMarkerColor(kOrange);
	g_off_ra30->SetMarkerColor(kGreen+2);
	g_amp_ra30->SetMarkerColor(kGreen+2);
	g_tau_ra30->SetMarkerColor(kGreen+2);
	g_off_ra_10_20_30->SetMarkerColor(kViolet);
	g_amp_ra_10_20_30->SetMarkerColor(kViolet);
	g_tau_ra_10_20_30->SetMarkerColor(kViolet);


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
		TFitResultPtr fit_kick_ra10 = kicks_ra10[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra10->SetPoint(g_off_ra10->GetN(),i+1,fit_kick_ra10->Parameter(2));
		g_off_ra10->SetPointError(g_off_ra10->GetN()-1,0,fit_kick_ra10->ParError(2));
		g_amp_ra10->SetPoint(g_amp_ra10->GetN(),i+1,fit_kick_ra10->Parameter(0));
		g_amp_ra10->SetPointError(g_amp_ra10->GetN()-1,0,fit_kick_ra10->ParError(0));
		g_tau_ra10->SetPoint(g_tau_ra10->GetN(),i+1,1000*fit_kick_ra10->Parameter(1));
		g_tau_ra10->SetPointError(g_tau_ra10->GetN()-1,0,1000*fit_kick_ra10->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra20 = kicks_ra20[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra20->SetPoint(g_off_ra20->GetN(),i+1,fit_kick_ra20->Parameter(2));
		g_off_ra20->SetPointError(g_off_ra20->GetN()-1,0,fit_kick_ra20->ParError(2));
		g_amp_ra20->SetPoint(g_amp_ra20->GetN(),i+1,fit_kick_ra20->Parameter(0));
		g_amp_ra20->SetPointError(g_amp_ra20->GetN()-1,0,fit_kick_ra20->ParError(0));
		g_tau_ra20->SetPoint(g_tau_ra20->GetN(),i+1,1000*fit_kick_ra20->Parameter(1));
		g_tau_ra20->SetPointError(g_tau_ra20->GetN()-1,0,1000*fit_kick_ra20->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra30 = kicks_ra30[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra30->SetPoint(g_off_ra30->GetN(),i+1,fit_kick_ra30->Parameter(2));
		g_off_ra30->SetPointError(g_off_ra30->GetN()-1,0,fit_kick_ra30->ParError(2));
		g_amp_ra30->SetPoint(g_amp_ra30->GetN(),i+1,fit_kick_ra30->Parameter(0));
		g_amp_ra30->SetPointError(g_amp_ra30->GetN()-1,0,fit_kick_ra30->ParError(0));
		g_tau_ra30->SetPoint(g_tau_ra30->GetN(),i+1,1000*fit_kick_ra30->Parameter(1));
		g_tau_ra30->SetPointError(g_tau_ra30->GetN()-1,0,1000*fit_kick_ra30->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra_10_20_30 = kicks_ra_10_20_30[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra_10_20_30->SetPoint(g_off_ra_10_20_30->GetN(),i+1,fit_kick_ra_10_20_30->Parameter(2));
		g_off_ra_10_20_30->SetPointError(g_off_ra_10_20_30->GetN()-1,0,fit_kick_ra_10_20_30->ParError(2));
		g_amp_ra_10_20_30->SetPoint(g_amp_ra_10_20_30->GetN(),i+1,fit_kick_ra_10_20_30->Parameter(0));
		g_amp_ra_10_20_30->SetPointError(g_amp_ra_10_20_30->GetN()-1,0,fit_kick_ra_10_20_30->ParError(0));
		g_tau_ra_10_20_30->SetPoint(g_tau_ra_10_20_30->GetN(),i+1,1000*fit_kick_ra_10_20_30->Parameter(1));
		g_tau_ra_10_20_30->SetPointError(g_tau_ra_10_20_30->GetN()-1,0,1000*fit_kick_ra_10_20_30->ParError(1));

		//kicks[i]->Add(f_exp,-1);
		//kicks_ra10[i]->Add(f_exp,-1);
		//kicks_ra20[i]->Add(f_exp,-1);
		//kicks_ra30[i]->Add(f_exp,-1);
		//kicks_ra_10_20_30[i]->Add(f_exp,-1);
	}

	TH1D** h1_fft_kicks = new TH1D* [8];
	TH1D** h1_fft_kicks_ra10 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra20 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra30 = new TH1D* [8];
	TH1D** h1_fft_kicks_ra_10_20_30 = new TH1D* [8];
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
		TH1 *fft_histogram_ra20 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra20 = SetupFFT(kicks_ra20[i], fft_xmin, fft_xmax);
		fft_histogram_ra20 = fftResidualInit_ra20->FFT(fft_histogram_ra20,"MAG");
		h1_fft_kicks_ra20[i] = RescaleAxis(fft_histogram_ra20, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra20[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra20[i]->SetStats(0);
		h1_fft_kicks_ra20[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra20[i]->Scale(1.0 / h1_fft_kicks_ra20[i]->Integral());
		h1_fft_kicks_ra20[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra20[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra30 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra30 = SetupFFT(kicks_ra30[i], fft_xmin, fft_xmax);
		fft_histogram_ra30 = fftResidualInit_ra30->FFT(fft_histogram_ra30,"MAG");
		h1_fft_kicks_ra30[i] = RescaleAxis(fft_histogram_ra30, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra30[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra30[i]->SetStats(0);
		h1_fft_kicks_ra30[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra30[i]->Scale(1.0 / h1_fft_kicks_ra30[i]->Integral());
		h1_fft_kicks_ra30[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra30[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra_10_20_30 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra_10_20_30 = SetupFFT(kicks_ra_10_20_30[i], fft_xmin, fft_xmax);
		fft_histogram_ra_10_20_30 = fftResidualInit_ra_10_20_30->FFT(fft_histogram_ra_10_20_30,"MAG");
		h1_fft_kicks_ra_10_20_30[i] = RescaleAxis(fft_histogram_ra_10_20_30, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra_10_20_30[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra_10_20_30[i]->SetStats(0);
		h1_fft_kicks_ra_10_20_30[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra_10_20_30[i]->Scale(1.0 / h1_fft_kicks_ra_10_20_30[i]->Integral());
		h1_fft_kicks_ra_10_20_30[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra_10_20_30[i]->GetYaxis()->SetRangeUser(0,0.01);
	}

	TCanvas* can = new TCanvas("can","",1400,800);
	can->Divide(5,2);
	can->cd(1);
	kicks[0]->Draw("HIST");
	can->cd(2);
	kicks_ra10[0]->Draw("HIST");
	can->cd(3);
	kicks_ra20[0]->Draw("HIST");
	can->cd(4);
	kicks_ra30[0]->Draw("HIST");
	can->cd(5);
	kicks_ra_10_20_30[0]->Draw("HIST");
	can->cd(6);
	h1_fft_kicks[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks: "<<h1_fft_kicks[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(7);
	h1_fft_kicks_ra10[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra10: "<<h1_fft_kicks_ra10[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(8);
	h1_fft_kicks_ra20[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra20: "<<h1_fft_kicks_ra20[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(9);
	h1_fft_kicks_ra30[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra30: "<<h1_fft_kicks_ra30[0]->Integral()<<"\n";
	gPad->SetLogx();
	can->cd(10);
	h1_fft_kicks_ra_10_20_30[0]->Draw("HIST");
	cout<<"FFT integral h1_fft_kicks_ra_10_20_30: "<<h1_fft_kicks_ra_10_20_30[0]->Integral()<<"\n";
	gPad->SetLogx();


	TCanvas* can8 = new TCanvas("can8","",1400,800);
	can8->Divide(5,2);
	can8->cd(1);
	kicks[7]->Draw("HIST");
	can8->cd(2);
	kicks_ra10[7]->Draw("HIST");
	can8->cd(3);
	kicks_ra20[7]->Draw("HIST");
	can8->cd(4);
	kicks_ra30[7]->Draw("HIST");
	can8->cd(5);
	kicks_ra_10_20_30[7]->Draw("HIST");
	can8->cd(6);
	h1_fft_kicks[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(7);
	h1_fft_kicks_ra10[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(8);
	h1_fft_kicks_ra20[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(9);
	h1_fft_kicks_ra30[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(10);
	h1_fft_kicks_ra_10_20_30[7]->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can_fits = new TCanvas("can_fits","",1500,500);
	can_fits->Divide(3,1);
	can_fits->cd(1);
	g_off->GetYaxis()->SetRangeUser(-1,1);
	g_off->Draw("APZ");
	g_off_ra10->Draw("PZ");
	g_off_ra20->Draw("PZ");
	g_off_ra30->Draw("PZ");
	g_off_ra_10_20_30->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(2);
	g_amp->GetYaxis()->SetRangeUser(0,25);
	g_amp->Draw("APZ");
	g_amp_ra10->Draw("PZ");
	g_amp_ra20->Draw("PZ");
	g_amp_ra30->Draw("PZ");
	g_amp_ra_10_20_30->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(3);
	g_tau->GetYaxis()->SetRangeUser(20,100);
	g_tau->Draw("APZ");
	g_tau_ra10->Draw("PZ");
	g_tau_ra20->Draw("PZ");
	g_tau_ra30->Draw("PZ");
	g_tau_ra_10_20_30->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);

}

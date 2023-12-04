#include "analysis_tools.C"


void test_smooth(TString input_file, TString output_file=""){

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


	TH1D* trace_smooth = (TH1D*)trace->Clone("trace_smooth");
	trace_smooth->Smooth();
	trace_smooth->SetTitle(Form("%s (Smooth)",trace_smooth->GetTitle()));
	TH1D** kicks_smooth = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_smooth[i] = (TH1D*)kicks[i]->Clone(Form("kick%d_smooth",i+1));
		kicks_smooth[i]->Smooth();
		kicks_smooth[i]->SetTitle(Form("%s (Smooth)",kicks_smooth[i]->GetTitle()));
	}
	TH1D* kick8long_smooth = (TH1D*)kick8long->Clone("kick8long_smooth");
	kick8long_smooth->Smooth();
	kick8long_smooth->SetTitle(Form("%s (Smooth)",kick8long_smooth->GetTitle()));

	int rebin=10;
	TH1D* trace_rebin = (TH1D*)trace->Clone("trace_rebin");
	trace_rebin->SetTitle(Form("%s (Rebin %d)",trace_rebin->GetTitle(),rebin));
	trace_rebin->Rebin(rebin);
	trace_rebin->Scale(1./rebin);
	TH1D** kicks_rebin = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_rebin[i] = (TH1D*)kicks[i]->Clone(Form("kick%d_rebin",i+1));
		kicks_rebin[i]->SetTitle(Form("%s (Rebin %d)",kicks_rebin[i]->GetTitle(),rebin));
		kicks_rebin[i]->Rebin(rebin);
		kicks_rebin[i]->Scale(1./rebin);
	}
	TH1D* kick8long_rebin = (TH1D*)kick8long->Clone("kick8long_rebin");
	kick8long_rebin->SetTitle(Form("%s (Rebin %d)",kick8long_rebin->GetTitle(),rebin));
	kick8long_rebin->Rebin(rebin);
	kick8long_rebin->Scale(1./rebin);


	double lowpass_freq = 50.0;
	TH1D* trace_lowpass = lowpass(trace,lowpass_freq);
	trace_lowpass->SetTitle(Form("%s (Lowpass %.1f kHz)",trace->GetTitle(),lowpass_freq));	
	TH1D** kicks_lowpass = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_lowpass[i] = lowpass(kicks[i],lowpass_freq);
		kicks_lowpass[i]->SetTitle(Form("%s (Lowpass %.1f kHz)",kicks[i]->GetTitle(),lowpass_freq));	
	}
	TH1D* kick8long_lowpass = lowpass(kick8long,lowpass_freq);
	kick8long_lowpass->SetTitle(Form("%s (Lowpass %.1f kHz)",kick8long->GetTitle(),lowpass_freq));	

	TH1D* trace_ra1020 = runningAverage(runningAverage(trace,10,true),20,true);
	trace_ra1020->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),1020));	
	TH1D** kicks_ra1020 = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra1020[i] = runningAverage(runningAverage(kicks[i],10,true),20,true);
		kicks_ra1020[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),1020));	
	}
	TH1D* kick8long_ra1020 = runningAverage(runningAverage(kick8long,10,true),20,true);
	kick8long_ra1020->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),1020));	


	new TCanvas();
	trace->SetLineWidth(2);
	trace_smooth->SetLineWidth(2);
	trace_rebin->SetLineWidth(2);
	trace_lowpass->SetLineWidth(2);
	trace_ra1020->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_smooth->SetLineColor(kRed);
	trace_rebin->SetLineColor(kOrange);
	trace_lowpass->SetLineColor(kGreen+2);
	trace_ra1020->SetLineColor(kViolet);
	trace->Draw("HIST");
	trace_smooth->Draw("HIST SAME");
	trace_rebin->Draw("HIST SAME");
	trace_lowpass->Draw("HIST SAME");
	trace_ra1020->Draw("HIST SAME");


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_smooth[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_rebin[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_lowpass[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra1020[i]->GetYaxis()->SetRangeUser(-60,60);
	}

	TCanvas* can_kicks = new TCanvas("can_kicks","",1600,800);
	can_kicks->Divide(2,1);
	can_kicks->cd(1);
	kicks[0]->SetLineWidth(2);
	kicks_smooth[0]->SetLineWidth(2);
	kicks_rebin[0]->SetLineWidth(2);
	kicks_lowpass[0]->SetLineWidth(2);
	kicks_ra1020[0]->SetLineWidth(2);
	kicks[0]->SetLineColor(kBlue);
	kicks_smooth[0]->SetLineColor(kRed);
	kicks_rebin[0]->SetLineColor(kOrange);
	kicks_lowpass[0]->SetLineColor(kGreen+2);
	kicks_ra1020[0]->SetLineColor(kViolet);
	kicks[0]->Draw("HIST");
	kicks_smooth[0]->Draw("HIST SAME");
	kicks_rebin[0]->Draw("HIST SAME");
	kicks_lowpass[0]->Draw("HIST SAME");
	kicks_ra1020[0]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks->cd(2);
	kicks[7]->SetLineWidth(2);
	kicks_smooth[7]->SetLineWidth(2);
	kicks_rebin[7]->SetLineWidth(2);
	kicks_lowpass[7]->SetLineWidth(2);
	kicks_ra1020[7]->SetLineWidth(2);
	kicks[7]->SetLineColor(kBlue);
	kicks_smooth[7]->SetLineColor(kRed);
	kicks_rebin[7]->SetLineColor(kOrange);
	kicks_lowpass[7]->SetLineColor(kGreen+2);
	kicks_ra1020[7]->SetLineColor(kViolet);
	kicks[7]->Draw("HIST");
	kicks_smooth[7]->Draw("HIST SAME");
	kicks_rebin[7]->Draw("HIST SAME");
	kicks_lowpass[7]->Draw("HIST SAME");
	kicks_ra1020[7]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);


	TH1D* kick1_zoom = (TH1D*)kicks[0]->Clone("kick1_zoom");
	TH1D* kick1_zoom_smooth = (TH1D*)kicks_smooth[0]->Clone("kick1_zoom_smooth");
	TH1D* kick1_zoom_rebin = (TH1D*)kicks_rebin[0]->Clone("kick1_zoom_rebin");
	TH1D* kick1_zoom_lowpass = (TH1D*)kicks_lowpass[0]->Clone("kick1_zoom_lowpass");
	TH1D* kick1_zoom_ra1020 = (TH1D*)kicks_ra1020[0]->Clone("kick1_zoom_ra1020");

	TH1D* kick8_zoom = (TH1D*)kicks[7]->Clone("kick8_zoom");
	TH1D* kick8_zoom_smooth = (TH1D*)kicks_smooth[7]->Clone("kick8_zoom_smooth");
	TH1D* kick8_zoom_rebin = (TH1D*)kicks_rebin[7]->Clone("kick8_zoom_rebin");
	TH1D* kick8_zoom_lowpass = (TH1D*)kicks_lowpass[7]->Clone("kick8_zoom_lowpass");
	TH1D* kick8_zoom_ra1020 = (TH1D*)kicks_ra1020[7]->Clone("kick8_zoom_ra1020");

	double trace_offset = 25;
	for (int bn=1; bn<=kick1_zoom->GetNbinsX(); bn++){
		kick1_zoom_smooth->SetBinContent(bn,kick1_zoom_smooth->GetBinContent(bn) - trace_offset);
		kick1_zoom_rebin->SetBinContent(bn,kick1_zoom_rebin->GetBinContent(bn) - 2*trace_offset);
		kick1_zoom_lowpass->SetBinContent(bn,kick1_zoom_lowpass->GetBinContent(bn) - 3*trace_offset);
		kick1_zoom_ra1020->SetBinContent(bn,kick1_zoom_ra1020->GetBinContent(bn) - 4*trace_offset);

		kick8_zoom_smooth->SetBinContent(bn,kick8_zoom_smooth->GetBinContent(bn) - trace_offset);
		kick8_zoom_rebin->SetBinContent(bn,kick8_zoom_rebin->GetBinContent(bn) - 2*trace_offset);
		kick8_zoom_lowpass->SetBinContent(bn,kick8_zoom_lowpass->GetBinContent(bn) - 3*trace_offset);
		kick8_zoom_ra1020->SetBinContent(bn,kick8_zoom_ra1020->GetBinContent(bn) - 4*trace_offset);
	}

	double x_zoom_from = 0.0;
	double x_zoom_to = 1.0;
	kick1_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_smooth->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_rebin->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_lowpass->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick1_zoom_ra1020->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_smooth->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_rebin->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_lowpass->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);
	kick8_zoom_ra1020->GetXaxis()->SetRangeUser(x_zoom_from,x_zoom_to);

	kick1_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_smooth->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_rebin->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_lowpass->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick1_zoom_ra1020->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_smooth->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_rebin->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_lowpass->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);
	kick8_zoom_ra1020->GetYaxis()->SetRangeUser(-5*trace_offset,trace_offset);

	TCanvas* can_kicks_zoom = new TCanvas("can_kicks_zoom","",1800,900);
	can_kicks_zoom->Divide(2,1);
	can_kicks_zoom->cd(1);
	kick1_zoom->Draw("HIST");
	kick1_zoom_smooth->Draw("HIST SAME");
	kick1_zoom_rebin->Draw("HIST SAME");
	kick1_zoom_lowpass->Draw("HIST SAME");
	kick1_zoom_ra1020->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks_zoom->cd(2);
	kick8_zoom->Draw("HIST");
	kick8_zoom_smooth->Draw("HIST SAME");
	kick8_zoom_rebin->Draw("HIST SAME");
	kick8_zoom_lowpass->Draw("HIST SAME");
	kick8_zoom_ra1020->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);



	//Fits
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.5,0.0);
	double blumlein_fit_start = -0.4;
	double blumlein_fit_end = -0.2; 

	f_blumlein->SetParameters(60.0,0.5*(blumlein_fit_start+blumlein_fit_end),-1000.0);
	TFitResultPtr fit1 = kicks[0]->Fit("f_blumlein","QNS","",blumlein_fit_start,blumlein_fit_end);
	double blum_1 = fit1->Parameter(0);

	f_blumlein->SetParameters(1.0,0.5*(blumlein_fit_start+blumlein_fit_end),-10.0);
	TFitResultPtr fit2 = kicks_smooth[0]->Fit("f_blumlein","QNS","",blumlein_fit_start,blumlein_fit_end);
	double blum_2 = fit2->Parameter(0);

	cout<<blum_1<<" "<<blum_2<<"\n";

	double ratio = blum_1/blum_2;

	trace_smooth->Scale(ratio);
	for (int i=0; i<8; i++){
		kicks_smooth[i]->Scale(ratio);
	}
	kick8long_smooth->Scale(ratio);

	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_smooth[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_rebin[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_lowpass[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra1020[i]->GetYaxis()->SetRangeUser(-60,60);
	}


	TGraphErrors* g_off = new TGraphErrors();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	TGraphErrors* g_off_smooth = new TGraphErrors();
	TGraphErrors* g_amp_smooth = new TGraphErrors();
	TGraphErrors* g_tau_smooth = new TGraphErrors();
	TGraphErrors* g_off_rebin = new TGraphErrors();
	TGraphErrors* g_amp_rebin = new TGraphErrors();
	TGraphErrors* g_tau_rebin = new TGraphErrors();
	TGraphErrors* g_off_lowpass = new TGraphErrors();
	TGraphErrors* g_amp_lowpass = new TGraphErrors();
	TGraphErrors* g_tau_lowpass = new TGraphErrors();
	TGraphErrors* g_off_ra1020 = new TGraphErrors();
	TGraphErrors* g_amp_ra1020 = new TGraphErrors();
	TGraphErrors* g_tau_ra1020 = new TGraphErrors();
	g_off->SetTitle("Offset;Kick;Offset [mV]");
	g_amp->SetTitle("Amplitude;Kick;Amplitude [mV]");
	g_tau->SetTitle("Lifetime;Kick;Lifetime [#mus]");
	g_off_smooth->SetTitle("Offset (Smooth);Kick;Offset [mV]");
	g_amp_smooth->SetTitle("Amplitude (Smooth);Kick;Amplitude [mV]");
	g_tau_smooth->SetTitle("Lifetime (Smooth);Kick;Lifetime [#mus]");
	g_off_rebin->SetTitle("Offset (Rebin);Kick;Offset [mV]");
	g_amp_rebin->SetTitle("Amplitude (Rebin);Kick;Amplitude [mV]");
	g_tau_rebin->SetTitle("Lifetime (Rebin);Kick;Lifetime [#mus]");
	g_off_lowpass->SetTitle("Offset (Lowpass);Kick;Offset [mV]");
	g_amp_lowpass->SetTitle("Amplitude (Lowpass);Kick;Amplitude [mV]");
	g_tau_lowpass->SetTitle("Lifetime (Lowpass);Kick;Lifetime [#mus]");
	g_off_ra1020->SetTitle("Offset (RA 1020);Kick;Offset [mV]");
	g_amp_ra1020->SetTitle("Amplitude (RA 1020);Kick;Amplitude [mV]");
	g_tau_ra1020->SetTitle("Lifetime (RA 1020);Kick;Lifetime [#mus]");
	g_off->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_off_smooth->SetMarkerStyle(20);
	g_amp_smooth->SetMarkerStyle(20);
	g_tau_smooth->SetMarkerStyle(20);
	g_off_rebin->SetMarkerStyle(20);
	g_amp_rebin->SetMarkerStyle(20);
	g_tau_rebin->SetMarkerStyle(20);
	g_off_lowpass->SetMarkerStyle(20);
	g_amp_lowpass->SetMarkerStyle(20);
	g_tau_lowpass->SetMarkerStyle(20);
	g_off_ra1020->SetMarkerStyle(20);
	g_amp_ra1020->SetMarkerStyle(20);
	g_tau_ra1020->SetMarkerStyle(20);
	g_off->SetMarkerColor(kBlue);
	g_amp->SetMarkerColor(kBlue);
	g_tau->SetMarkerColor(kBlue);
	g_off_smooth->SetMarkerColor(kRed);
	g_amp_smooth->SetMarkerColor(kRed);
	g_tau_smooth->SetMarkerColor(kRed);
	g_off_rebin->SetMarkerColor(kOrange);
	g_amp_rebin->SetMarkerColor(kOrange);
	g_tau_rebin->SetMarkerColor(kOrange);
	g_off_lowpass->SetMarkerColor(kGreen+2);
	g_amp_lowpass->SetMarkerColor(kGreen+2);
	g_tau_lowpass->SetMarkerColor(kGreen+2);
	g_off_ra1020->SetMarkerColor(kViolet);
	g_amp_ra1020->SetMarkerColor(kViolet);
	g_tau_ra1020->SetMarkerColor(kViolet);


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
		TFitResultPtr fit_kick_smooth = kicks_smooth[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_smooth->SetPoint(g_off_smooth->GetN(),i+1,fit_kick_smooth->Parameter(2));
		g_off_smooth->SetPointError(g_off_smooth->GetN()-1,0,fit_kick_smooth->ParError(2));
		g_amp_smooth->SetPoint(g_amp_smooth->GetN(),i+1,fit_kick_smooth->Parameter(0));
		g_amp_smooth->SetPointError(g_amp_smooth->GetN()-1,0,fit_kick_smooth->ParError(0));
		g_tau_smooth->SetPoint(g_tau_smooth->GetN(),i+1,1000*fit_kick_smooth->Parameter(1));
		g_tau_smooth->SetPointError(g_tau_smooth->GetN()-1,0,1000*fit_kick_smooth->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_rebin = kicks_rebin[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_rebin->SetPoint(g_off_rebin->GetN(),i+1,fit_kick_rebin->Parameter(2));
		g_off_rebin->SetPointError(g_off_rebin->GetN()-1,0,fit_kick_rebin->ParError(2));
		g_amp_rebin->SetPoint(g_amp_rebin->GetN(),i+1,fit_kick_rebin->Parameter(0));
		g_amp_rebin->SetPointError(g_amp_rebin->GetN()-1,0,fit_kick_rebin->ParError(0));
		g_tau_rebin->SetPoint(g_tau_rebin->GetN(),i+1,1000*fit_kick_rebin->Parameter(1));
		g_tau_rebin->SetPointError(g_tau_rebin->GetN()-1,0,1000*fit_kick_rebin->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_lowpass = kicks_lowpass[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_lowpass->SetPoint(g_off_lowpass->GetN(),i+1,fit_kick_lowpass->Parameter(2));
		g_off_lowpass->SetPointError(g_off_lowpass->GetN()-1,0,fit_kick_lowpass->ParError(2));
		g_amp_lowpass->SetPoint(g_amp_lowpass->GetN(),i+1,fit_kick_lowpass->Parameter(0));
		g_amp_lowpass->SetPointError(g_amp_lowpass->GetN()-1,0,fit_kick_lowpass->ParError(0));
		g_tau_lowpass->SetPoint(g_tau_lowpass->GetN(),i+1,1000*fit_kick_lowpass->Parameter(1));
		g_tau_lowpass->SetPointError(g_tau_lowpass->GetN()-1,0,1000*fit_kick_lowpass->ParError(1));

		f_exp->SetParameters(10,60,0);
		TFitResultPtr fit_kick_ra1020 = kicks_ra1020[i]->Fit(f_exp,"NS","",0.02,0.7);
		g_off_ra1020->SetPoint(g_off_ra1020->GetN(),i+1,fit_kick_ra1020->Parameter(2));
		g_off_ra1020->SetPointError(g_off_ra1020->GetN()-1,0,fit_kick_ra1020->ParError(2));
		g_amp_ra1020->SetPoint(g_amp_ra1020->GetN(),i+1,fit_kick_ra1020->Parameter(0));
		g_amp_ra1020->SetPointError(g_amp_ra1020->GetN()-1,0,fit_kick_ra1020->ParError(0));
		g_tau_ra1020->SetPoint(g_tau_ra1020->GetN(),i+1,1000*fit_kick_ra1020->Parameter(1));
		g_tau_ra1020->SetPointError(g_tau_ra1020->GetN()-1,0,1000*fit_kick_ra1020->ParError(1));

		//kicks[i]->Add(f_exp,-1);
		//kicks_smooth[i]->Add(f_exp,-1);
		//kicks_rebin[i]->Add(f_exp,-1);
		//kicks_lowpass[i]->Add(f_exp,-1);
		//kicks_ra1020[i]->Add(f_exp,-1);
	}

	TH1D** h1_fft_kicks = new TH1D* [8];
	TH1D** h1_fft_kicks_smooth = new TH1D* [8];
	TH1D** h1_fft_kicks_rebin = new TH1D* [8];
	TH1D** h1_fft_kicks_lowpass = new TH1D* [8];
	TH1D** h1_fft_kicks_ra1020 = new TH1D* [8];
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
		TH1 *fft_histogram_smooth = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_smooth = SetupFFT(kicks_smooth[i], fft_xmin, fft_xmax);
		fft_histogram_smooth = fftResidualInit_smooth->FFT(fft_histogram_smooth,"MAG");
		h1_fft_kicks_smooth[i] = RescaleAxis(fft_histogram_smooth, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_smooth[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_smooth[i]->SetStats(0);
		h1_fft_kicks_smooth[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_smooth[i]->Scale(1.0 / h1_fft_kicks_smooth[i]->Integral());
		h1_fft_kicks_smooth[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_smooth[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_rebin = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_rebin = SetupFFT(kicks_rebin[i], fft_xmin, fft_xmax);
		fft_histogram_rebin = fftResidualInit_rebin->FFT(fft_histogram_rebin,"MAG");
		h1_fft_kicks_rebin[i] = RescaleAxis(fft_histogram_rebin, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_rebin[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_rebin[i]->SetStats(0);
		h1_fft_kicks_rebin[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_rebin[i]->Scale(1.0 / h1_fft_kicks_rebin[i]->Integral());
		h1_fft_kicks_rebin[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),(0.5/dt)/rebin);
		h1_fft_kicks_rebin[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_lowpass = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_lowpass = SetupFFT(kicks_lowpass[i], fft_xmin, fft_xmax);
		fft_histogram_lowpass = fftResidualInit_lowpass->FFT(fft_histogram_lowpass,"MAG");
		h1_fft_kicks_lowpass[i] = RescaleAxis(fft_histogram_lowpass, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_lowpass[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_lowpass[i]->SetStats(0);
		h1_fft_kicks_lowpass[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_lowpass[i]->Scale(1.0 / h1_fft_kicks_lowpass[i]->Integral());
		h1_fft_kicks_lowpass[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_lowpass[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		TH1 *fft_histogram_ra1020 = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra1020 = SetupFFT(kicks_ra1020[i], fft_xmin, fft_xmax);
		fft_histogram_ra1020 = fftResidualInit_ra1020->FFT(fft_histogram_ra1020,"MAG");
		h1_fft_kicks_ra1020[i] = RescaleAxis(fft_histogram_ra1020, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra1020[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra1020[i]->SetStats(0);
		h1_fft_kicks_ra1020[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra1020[i]->Scale(1.0 / h1_fft_kicks_ra1020[i]->Integral());
		h1_fft_kicks_ra1020[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra1020[i]->GetYaxis()->SetRangeUser(0,0.01);
	}

	TCanvas* can = new TCanvas("can","",1400,800);
	can->Divide(5,2);
	can->cd(1);
	kicks[0]->Draw("HIST");
	can->cd(2);
	kicks_smooth[0]->Draw("HIST");
	can->cd(3);
	kicks_rebin[0]->Draw("HIST");
	can->cd(4);
	kicks_lowpass[0]->Draw("HIST");
	can->cd(5);
	kicks_ra1020[0]->Draw("HIST");
	can->cd(6);
	h1_fft_kicks[0]->Draw("HIST");
	gPad->SetLogx();
	can->cd(7);
	h1_fft_kicks_smooth[0]->Draw("HIST");
	gPad->SetLogx();
	can->cd(8);
	h1_fft_kicks_rebin[0]->Draw("HIST");
	gPad->SetLogx();
	can->cd(9);
	h1_fft_kicks_lowpass[0]->Draw("HIST");
	gPad->SetLogx();
	can->cd(10);
	h1_fft_kicks_ra1020[0]->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can8 = new TCanvas("can8","",1400,800);
	can8->Divide(5,2);
	can8->cd(1);
	kicks[7]->Draw("HIST");
	can8->cd(2);
	kicks_smooth[7]->Draw("HIST");
	can8->cd(3);
	kicks_rebin[7]->Draw("HIST");
	can8->cd(4);
	kicks_lowpass[7]->Draw("HIST");
	can8->cd(5);
	kicks_ra1020[7]->Draw("HIST");
	can8->cd(6);
	h1_fft_kicks[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(7);
	h1_fft_kicks_smooth[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(8);
	h1_fft_kicks_rebin[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(9);
	h1_fft_kicks_lowpass[7]->Draw("HIST");
	gPad->SetLogx();
	can8->cd(10);
	h1_fft_kicks_ra1020[7]->Draw("HIST");
	gPad->SetLogx();

	TCanvas* can_fits = new TCanvas("can_fits","",1500,500);
	can_fits->Divide(3,1);
	can_fits->cd(1);
	g_off->GetYaxis()->SetRangeUser(-1,1);
	g_off->Draw("APZ");
	g_off_smooth->Draw("PZ");
	g_off_rebin->Draw("PZ");
	g_off_lowpass->Draw("PZ");
	g_off_ra1020->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(2);
	g_amp->GetYaxis()->SetRangeUser(0,25);
	g_amp->Draw("APZ");
	g_amp_smooth->Draw("PZ");
	g_amp_rebin->Draw("PZ");
	g_amp_lowpass->Draw("PZ");
	g_amp_ra1020->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);
	can_fits->cd(3);
	g_tau->GetYaxis()->SetRangeUser(20,100);
	g_tau->Draw("APZ");
	g_tau_smooth->Draw("PZ");
	g_tau_rebin->Draw("PZ");
	g_tau_lowpass->Draw("PZ");
	g_tau_ra1020->Draw("PZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.5,0.7,0.9,0.9);

}

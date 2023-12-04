#include "analysis_tools.C"

void fit_eddycurrents(TString input_file, TString output_file=""){

	TFile* f = TFile::Open(input_file);


	bool saveOutput = false;
	if (output_file != ""){
		saveOutput = true;
	}


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


	TH1D* trace_ra = runningAverage(runningAverage(runningAverage(trace,5,true),10,true),15,true);
	trace_ra->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),51015));	
	TH1D** kicks_ra = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra[i] = runningAverage(runningAverage(runningAverage(kicks[i],5,true),10,true),15,true);
		kicks_ra[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),51015));	
	}
	TH1D* kick8long_ra = runningAverage(runningAverage(runningAverage(kick8long,5,true),10,true),15,true);
	kick8long_ra->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),51015));	



	new TCanvas();
	trace->SetLineWidth(2);
	trace_ra->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_ra->SetLineColor(kRed);
	trace->Draw("HIST");
	trace_ra->Draw("HIST SAME");


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra[i]->GetYaxis()->SetRangeUser(-60,60);
	}


	TCanvas* can_kicks = new TCanvas("can_kicks","",1600,800);
	can_kicks->Divide(2,1);
	can_kicks->cd(1);
	kicks[0]->SetLineWidth(2);
	kicks_ra[0]->SetLineWidth(2);
	kicks[0]->SetLineColor(kBlue);
	kicks_ra[0]->SetLineColor(kRed);
	kicks[0]->Draw("HIST");
	kicks_ra[0]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks->cd(2);
	kicks[7]->SetLineWidth(2);
	kicks_ra[7]->SetLineWidth(2);
	kicks[7]->SetLineColor(kBlue);
	kicks_ra[7]->SetLineColor(kRed);
	kicks[7]->Draw("HIST");
	kicks_ra[7]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);



	TH1D** h1_fft_kicks = new TH1D* [8];
	TH1D** h1_fft_kicks_ra = new TH1D* [8];
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
		TH1 *fft_histogram_ra = 0;
		TVirtualFFT::SetTransform(0);
		TH1D* fftResidualInit_ra = SetupFFT(kicks_ra[i], fft_xmin, fft_xmax);
		fft_histogram_ra = fftResidualInit_ra->FFT(fft_histogram_ra,"MAG");
		h1_fft_kicks_ra[i] = RescaleAxis(fft_histogram_ra, 1./(fft_xmax - fft_xmin));
		h1_fft_kicks_ra[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
		h1_fft_kicks_ra[i]->SetStats(0);
		h1_fft_kicks_ra[i]->SetName(Form("residualFFT_%d",i));
		h1_fft_kicks_ra[i]->Scale(1.0 / h1_fft_kicks_ra[i]->Integral());
		h1_fft_kicks_ra[i]->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
		h1_fft_kicks_ra[i]->GetYaxis()->SetRangeUser(0,0.01);
	}


	TText* text = new TText();
	TCanvas* can = new TCanvas("can","",1400,800);
	can->Divide(2,2);
	can->cd(1);
	kicks[0]->Draw("HIST");
	can->cd(2);
	kicks_ra[0]->Draw("HIST");
	can->cd(3);
	h1_fft_kicks[0]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks[0]->Integral()));
	gPad->SetLogx();
	can->cd(4);
	h1_fft_kicks_ra[0]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks_ra[0]->Integral()));
	gPad->SetLogx();


	TCanvas* can8 = new TCanvas("can8","",1400,800);
	can8->Divide(2,2);
	can8->cd(1);
	kicks[7]->Draw("HIST");
	can8->cd(2);
	kicks_ra[7]->Draw("HIST");
	can8->cd(3);
	h1_fft_kicks[7]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks[7]->Integral()));
	gPad->SetLogx();
	can8->cd(4);
	h1_fft_kicks_ra[7]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks_ra[7]->Integral()));
	gPad->SetLogx();






	//Fits

	TH1D** kicks_fit_exp = new TH1D* [8];
	TH1D** kicks_fit_exp_res = new TH1D* [8];
	for (int i=0; i<8; i++){
		kicks_fit_exp[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp",i+1));
		kicks_fit_exp_res[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp_res",i+1));
		
		kicks_fit_exp[i]->SetLineColor(kBlue);
		kicks_fit_exp_res[i]->SetLineColor(kBlue);
		kicks_fit_exp[i]->SetLineWidth(1);
		kicks_fit_exp_res[i]->SetLineWidth(1);
		kicks_fit_exp[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp_res[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp[i]->GetYaxis()->SetRangeUser(-20,20);
		kicks_fit_exp_res[i]->GetYaxis()->SetRangeUser(-20,20);
	}


	TCanvas* can_fits1 = new TCanvas("can_fits1","",800,800);
	can_fits1->Divide(2,4);
	TCanvas* can_fits1_res = new TCanvas("can_fits1_res","",800,800);
	can_fits1_res->Divide(2,4);


	TGraphErrors* g_off = new TGraphErrors();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	g_off->SetTitle("Offset;Kick;Offset [mV]");
	g_amp->SetTitle("Amplitude;Kick;Amplitude [mV]");
	g_tau->SetTitle("Lifetime;Kick;Lifetime [#mus]");
	g_off->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_off->SetMarkerColor(kRed);
	g_amp->SetMarkerColor(kRed);
	g_tau->SetMarkerColor(kRed);

	TF1* f_exp = new TF1("f_exp","[2]-[0]*exp(-x/[1])",0.02,0.7);
	f_exp->SetParameters(10,0.07,0);

	TF1* f_exp_17 = new TF1("f_exp_17","[2]-[0]*exp(-x/[1])*(",0.02,0.7);
	f_exp_17->SetParameters(10,0.07,0);

	for (int i=0; i<8; i++){
		can_fits1->cd(i+1);
		kicks_fit_exp[i]->Draw();
		kicks_fit_exp[i]->GetYaxis()->SetRangeUser(-20,20);

		cout<<"Fitting kick "<<i+1<<"\n";
		f_exp->SetParameters(10,0.07,0);
		f_exp->FixParameter(2,0.0);
		TFitResultPtr fit_kick = kicks_fit_exp[i]->Fit(f_exp,"S","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick->ParError(1));

		kicks_fit_exp_res[i]->Add(f_exp,-1);
		can_fits1_res->cd(i+1);
		kicks_fit_exp_res[i]->Draw();
		kicks_fit_exp_res[i]->GetYaxis()->SetRangeUser(-20,20);
	}


	TCanvas* can_fit_res1 = new TCanvas("can_fit_res1","",1500,500);
	can_fit_res1->Divide(3,1);
	can_fit_res1->cd(1);
	g_off->GetYaxis()->SetRangeUser(-1,1);
	g_off->Draw("APZ");
	gPad->SetGridy();
	can_fit_res1->cd(2);
	g_amp->GetYaxis()->SetRangeUser(0,25);
	g_amp->Draw("APZ");
	gPad->SetGridy();
	can_fit_res1->cd(3);
	g_tau->GetYaxis()->SetRangeUser(20,100);
	g_tau->Draw("APZ");
	gPad->SetGridy();

}
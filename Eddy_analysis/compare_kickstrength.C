#include "../analysis_tools.C"



void compare_kickstrength(){

	TFile* f1 = TFile::Open("fits/fitted_normalized_analysis_SD_R0_eddy_oct26_H25_B100_k194.root");
	TFile* f2 = TFile::Open("fits/fitted_normalized_analysis_SD_R0_eddy_oct25_H25_B100_k389.root");
	TFile* f3 = TFile::Open("fits/fitted_normalized_analysis_SD_R0_eddy_oct25_H25_B100_k583.root");
	TFile* f4 = TFile::Open("fits/fitted_normalized_analysis_SD_R0_eddy_oct24_H25_B100_k777.root");
	vector<double> strengths = {1.94, 3.89, 5.83, 7.77};

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;

	TH1D* f1_kick1 = (TH1D*)f1->Get("trace_kick1_ABnormalized_px");
	TH1D* f2_kick1 = (TH1D*)f2->Get("trace_kick1_ABnormalized_px");
	TH1D* f3_kick1 = (TH1D*)f3->Get("trace_kick1_ABnormalized_px");
	TH1D* f4_kick1 = (TH1D*)f4->Get("trace_kick1_ABnormalized_px");

	TH1D* f1_kick1_ra = smoothing(f1_kick1,"f1_kick1_ra");
	TH1D* f2_kick1_ra = smoothing(f2_kick1,"f2_kick1_ra");
	TH1D* f3_kick1_ra = smoothing(f3_kick1,"f3_kick1_ra");
	TH1D* f4_kick1_ra = smoothing(f4_kick1,"f4_kick1_ra");

	TH1D* f1_kick8 = (TH1D*)f1->Get("trace_kick8_ABnormalized_px");
	TH1D* f2_kick8 = (TH1D*)f2->Get("trace_kick8_ABnormalized_px");
	TH1D* f3_kick8 = (TH1D*)f3->Get("trace_kick8_ABnormalized_px");
	TH1D* f4_kick8 = (TH1D*)f4->Get("trace_kick8_ABnormalized_px");

	TH1D* f1_kick1_exp = (TH1D*)f1->Get("kick1_fit_exp");
	TH1D* f2_kick1_exp = (TH1D*)f2->Get("kick1_fit_exp");
	TH1D* f3_kick1_exp = (TH1D*)f3->Get("kick1_fit_exp");
	TH1D* f4_kick1_exp = (TH1D*)f4->Get("kick1_fit_exp");

	TH1D* f1_kick1_exp_res = (TH1D*)f1->Get("kick1_fit_exp_res");
	TH1D* f2_kick1_exp_res = (TH1D*)f2->Get("kick1_fit_exp_res");
	TH1D* f3_kick1_exp_res = (TH1D*)f3->Get("kick1_fit_exp_res");
	TH1D* f4_kick1_exp_res = (TH1D*)f4->Get("kick1_fit_exp_res");

	TH1D* f1_kick8_exp = (TH1D*)f1->Get("kick8_fit_exp");
	TH1D* f2_kick8_exp = (TH1D*)f2->Get("kick8_fit_exp");
	TH1D* f3_kick8_exp = (TH1D*)f3->Get("kick8_fit_exp");
	TH1D* f4_kick8_exp = (TH1D*)f4->Get("kick8_fit_exp");

	TH1D* f1_kick1_FFT = (TH1D*)f1->Get("trace_kick1_ABnormalized_px_FFT");
	TH1D* f2_kick1_FFT = (TH1D*)f2->Get("trace_kick1_ABnormalized_px_FFT");
	TH1D* f3_kick1_FFT = (TH1D*)f3->Get("trace_kick1_ABnormalized_px_FFT");
	TH1D* f4_kick1_FFT = (TH1D*)f4->Get("trace_kick1_ABnormalized_px_FFT");


	f1_kick1->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[0]));
	f2_kick1->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[1]));
	f3_kick1->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[2]));
	f4_kick1->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[3]));

	f1_kick1_ra->SetTitle(Form("Kick %.f %%",100.*strengths[0]/7.77));
	f2_kick1_ra->SetTitle(Form("Kick %.f %%",100.*strengths[1]/7.77));
	f3_kick1_ra->SetTitle(Form("Kick %.f %%",100.*strengths[2]/7.77));
	f4_kick1_ra->SetTitle(Form("Kick %.f %%",100.*strengths[3]/7.77));

	f1_kick8->SetTitle(Form("Kick 8 (Kick %.2f kV)",strengths[0]));
	f2_kick8->SetTitle(Form("Kick 8 (Kick %.2f kV)",strengths[1]));
	f3_kick8->SetTitle(Form("Kick 8 (Kick %.2f kV)",strengths[2]));
	f4_kick8->SetTitle(Form("Kick 8 (Kick %.2f kV)",strengths[3]));

	f1_kick1_FFT->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[0]));
	f2_kick1_FFT->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[1]));
	f3_kick1_FFT->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[2]));
	f4_kick1_FFT->SetTitle(Form("Kick 1 (Kick %.2f kV)",strengths[3]));




	TF1* f_vibration1 = new TF1("f_vibration1","[0]+[1]*sin([2]*x+[3])");
	TF1* f_vibration2 = new TF1("f_vibration2","[0]+[1]*sin([2]*x+[3])");
	TF1* f_vibration3 = new TF1("f_vibration3","[0]+[1]*sin([2]*x+[3])");
	TF1* f_vibration4 = new TF1("f_vibration4","[0]+[1]*sin([2]*x+[3])");

	f_vibration1->SetParameters(0,20,109,0);
	f_vibration1->SetParLimits(1,0,50);
	new TCanvas();
	f1_kick1_ra->Draw("HIST");
	TFitResultPtr fit_vibration1 = f1_kick1_ra->Fit("f_vibration1","QS+","",2,4);
	f1_kick1_ra->GetFunction("f_vibration1")->Draw("SAME");
	f_vibration2->SetParameters(0,20,109,0);
	f_vibration2->SetParLimits(1,0,50);
	new TCanvas();
	f2_kick1_ra->Draw("HIST");
	TFitResultPtr fit_vibration2 = f2_kick1_ra->Fit("f_vibration2","QS+","",2,4);
	f2_kick1_ra->GetFunction("f_vibration2")->Draw("SAME");
	f_vibration3->SetParameters(0,20,109,0);
	f_vibration3->SetParLimits(1,0,50);
	new TCanvas();
	f3_kick1_ra->Draw("HIST");
	TFitResultPtr fit_vibration3 = f3_kick1_ra->Fit("f_vibration3","QS+","",2,4);
	f3_kick1_ra->GetFunction("f_vibration3")->Draw("SAME");
	f_vibration4->SetParameters(0,20,109,0);
	f_vibration4->SetParLimits(1,0,50);
	new TCanvas();
	f4_kick1_ra->Draw("HIST");
	TFitResultPtr fit_vibration4 = f4_kick1_ra->Fit("f_vibration4","QS+","",2,4);
	f4_kick1_ra->GetFunction("f_vibration4")->Draw("SAME");

	TGraphErrors* g_vibration = new TGraphErrors();
	if (fit_vibration1 >= 0){
		g_vibration->SetPoint(0,100.*strengths[0]/7.77,fit_vibration1->Parameter(1));
		g_vibration->SetPointError(0,0,fit_vibration1->ParError(1));
	}
	if (fit_vibration2 >= 0){
		g_vibration->SetPoint(1,100.*strengths[1]/7.77,fit_vibration2->Parameter(1));
		g_vibration->SetPointError(1,0,fit_vibration2->ParError(1));
	}
	if (fit_vibration3 >= 0){
		g_vibration->SetPoint(2,100.*strengths[2]/7.77,fit_vibration3->Parameter(1));
		g_vibration->SetPointError(2,0,fit_vibration3->ParError(1));
	}
	if (fit_vibration4 >= 0){
		g_vibration->SetPoint(3,100.*strengths[3]/7.77,fit_vibration4->Parameter(1));
		g_vibration->SetPointError(3,0,fit_vibration4->ParError(1));
	}
	//new TCanvas();
	g_vibration->Sort();
	g_vibration->SetName("g_vibration");
	g_vibration->SetTitle("Vibration amplitude");
	g_vibration->GetXaxis()->SetTitle("Kicker strength [%]");
	g_vibration->GetYaxis()->SetTitle("Amplitude [mV]");
	g_vibration->SetMarkerStyle(20);
	//g_vibration->Draw("APL");


	TGraphErrors* g_blum_amp = new TGraphErrors();
	double amp1 = f1_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr1 = f1_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base1 = f1_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr1 = f1_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak1 = f1_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr1 = f1_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum1 = peak1 - base1;
	double blumerr1 = sqrt(baseerr1*baseerr1 + peakerr1*peakerr1);
	g_blum_amp->SetPoint(0,blum1,amp1);
	g_blum_amp->SetPointError(0,blumerr1,amperr1);
	double amp2 = f2_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr2 = f2_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base2 = f2_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr2 = f2_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak2 = f2_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr2 = f2_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum2 = peak2 - base2;
	double blumerr2 = sqrt(baseerr2*baseerr2 + peakerr2*peakerr2);
	g_blum_amp->SetPoint(1,blum2,amp2);
	g_blum_amp->SetPointError(1,blumerr2,amperr2);
	double amp3 = f3_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr3 = f3_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base3 = f3_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr3 = f3_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak3 = f3_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr3 = f3_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum3 = peak3 - base3;
	double blumerr3 = sqrt(baseerr3*baseerr3 + peakerr3*peakerr3);
	g_blum_amp->SetPoint(2,blum3,amp3);
	g_blum_amp->SetPointError(2,blumerr3,amperr3);
	double amp4 = f4_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr4 = f4_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base4 = f4_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr4 = f4_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak4 = f4_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr4 = f4_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum4 = peak4 - base4;
	double blumerr4 = sqrt(baseerr4*baseerr4 + peakerr4*peakerr4);
	g_blum_amp->SetPoint(3,blum4,amp4);
	g_blum_amp->SetPointError(3,blumerr4,amperr4);

	//new TCanvas();
	g_blum_amp->GetXaxis()->SetTitle("Blumlein [mV]");
	g_blum_amp->GetYaxis()->SetTitle("Exp amplitude [mV]");
	g_blum_amp->SetMarkerStyle(20);
	//g_blum_amp->Draw("APL");



	TGraphErrors* g_kick_blum = new TGraphErrors();
	g_kick_blum->SetPoint(g_kick_blum->GetN(),100.*strengths[0]/7.77,blum1);
	g_kick_blum->SetPointError(g_kick_blum->GetN()-1,0,blumerr1);
	g_kick_blum->SetPoint(g_kick_blum->GetN(),100.*strengths[1]/7.77,blum2);
	g_kick_blum->SetPointError(g_kick_blum->GetN()-1,0,blumerr2);
	g_kick_blum->SetPoint(g_kick_blum->GetN(),100.*strengths[2]/7.77,blum3);
	g_kick_blum->SetPointError(g_kick_blum->GetN()-1,0,blumerr3);
	g_kick_blum->SetPoint(g_kick_blum->GetN(),100.*strengths[3]/7.77,blum4);
	g_kick_blum->SetPointError(g_kick_blum->GetN()-1,0,blumerr4);
	//new TCanvas();
	g_kick_blum->SetTitle("Blumlein amplitude");
	g_kick_blum->GetXaxis()->SetTitle("Kicker strength [%]");
	g_kick_blum->GetYaxis()->SetTitle("Amplitude [mV]");
	g_kick_blum->SetMarkerStyle(20);
	//g_kick_blum->Draw("APL");



	TGraphErrors* g_kick_amp = new TGraphErrors();
	g_kick_amp->SetPoint(g_kick_amp->GetN(),100.*strengths[0]/7.77,amp1);
	g_kick_amp->SetPointError(g_kick_amp->GetN()-1,0,amperr1);
	g_kick_amp->SetPoint(g_kick_amp->GetN(),100.*strengths[1]/7.77,amp2);
	g_kick_amp->SetPointError(g_kick_amp->GetN()-1,0,amperr2);
	g_kick_amp->SetPoint(g_kick_amp->GetN(),100.*strengths[2]/7.77,amp3);
	g_kick_amp->SetPointError(g_kick_amp->GetN()-1,0,amperr3);
	g_kick_amp->SetPoint(g_kick_amp->GetN(),100.*strengths[3]/7.77,amp4);
	g_kick_amp->SetPointError(g_kick_amp->GetN()-1,0,amperr4);
	//new TCanvas();
	g_kick_amp->SetTitle("Transient amplitude");
	g_kick_amp->GetXaxis()->SetTitle("Kicker strength [%]");
	g_kick_amp->GetYaxis()->SetTitle("Amplitude [mV]");
	g_kick_amp->SetMarkerStyle(20);
	//g_kick_amp->Draw("APL");


	//normalize
	/*
	double norm_blum = g_kick_blum->GetPointY(g_kick_blum->GetN()-1);
	for (int i=0; i<g_kick_blum->GetN(); i++){
		g_kick_blum->SetPointY(i,g_kick_blum->GetPointY(i)/norm_blum);
	}
	double norm_amp = g_kick_amp->GetPointY(g_kick_amp->GetN()-1);
	for (int i=0; i<g_kick_amp->GetN(); i++){
		g_kick_amp->SetPointY(i,g_kick_amp->GetPointY(i)/norm_amp);
	}
	double norm_vibration = g_vibration->GetPointY(g_vibration->GetN()-1);
	for (int i=0; i<g_vibration->GetN(); i++){
		g_vibration->SetPointY(i,g_vibration->GetPointY(i)/norm_vibration);
	}
	*/
	new TCanvas();
	g_kick_blum->SetMarkerColor(kBlack);
	g_kick_blum->SetLineColor(kBlack);
	g_kick_amp->SetMarkerColor(kBlue);
	g_kick_amp->SetLineColor(kBlue);
	g_vibration->SetMarkerColor(kRed);
	g_vibration->SetLineColor(kRed);
	g_kick_blum->GetXaxis()->SetRangeUser(0,105);
	g_kick_blum->GetYaxis()->SetRangeUser(0,60);
	g_kick_blum->Draw("APL");
	g_kick_amp->Draw("PL");
	g_vibration->Draw("PL");
	gPad->BuildLegend();


	int bin17 = f1_kick1_FFT->FindBin(17);
	int bin18 = f1_kick1_FFT->FindBin(18);
	TGraphErrors* g_amp_17kHz = new TGraphErrors();
	double fft1_17_err;
	double fft1_17 = f1_kick1_FFT->IntegralAndError(bin17,bin18,fft1_17_err);
	g_amp_17kHz->SetPoint(0,amp1,fft1_17);
	g_amp_17kHz->SetPointError(0,amperr1,fft1_17_err);
	double fft2_17_err;
	double fft2_17 = f2_kick1_FFT->IntegralAndError(bin17,bin18,fft2_17_err);
	g_amp_17kHz->SetPoint(1,amp2,fft2_17);
	g_amp_17kHz->SetPointError(1,amperr2,fft2_17_err);
	double fft3_17_err;
	double fft3_17 = f3_kick1_FFT->IntegralAndError(bin17,bin18,fft3_17_err);
	g_amp_17kHz->SetPoint(2,amp3,fft3_17);
	g_amp_17kHz->SetPointError(2,amperr3,fft3_17_err);
	double fft4_17_err;
	double fft4_17 = f4_kick1_FFT->IntegralAndError(bin17,bin18,fft4_17_err);
	g_amp_17kHz->SetPoint(3,amp4,fft4_17);
	g_amp_17kHz->SetPointError(3,amperr4,fft4_17_err);

	new TCanvas();
	g_amp_17kHz->GetXaxis()->SetTitle("Exp amplitude [mV]");
	g_amp_17kHz->GetYaxis()->SetTitle("FFT 17kHz amplitude [arb. u.]");
	g_amp_17kHz->SetMarkerStyle(20);
	g_amp_17kHz->Draw("APL");


	int bin4 = f1_kick1_FFT->FindBin(4);
	int bin4p5 = f1_kick1_FFT->FindBin(4.5);
	TGraphErrors* g_amp_4kHz = new TGraphErrors();
	double fft1_4_err;
	double fft1_4 = f1_kick1_FFT->IntegralAndError(bin4,bin18,fft1_4_err);
	g_amp_4kHz->SetPoint(0,amp1,fft1_4);
	g_amp_4kHz->SetPointError(0,amperr1,fft1_4_err);
	double fft2_4_err;
	double fft2_4 = f2_kick1_FFT->IntegralAndError(bin4,bin18,fft2_4_err);
	g_amp_4kHz->SetPoint(1,amp2,fft2_4);
	g_amp_4kHz->SetPointError(1,amperr2,fft2_4_err);
	double fft3_4_err;
	double fft3_4 = f3_kick1_FFT->IntegralAndError(bin4,bin18,fft3_4_err);
	g_amp_4kHz->SetPoint(2,amp3,fft3_4);
	g_amp_4kHz->SetPointError(2,amperr3,fft3_4_err);
	double fft4_4_err;
	double fft4_4 = f4_kick1_FFT->IntegralAndError(bin4,bin18,fft4_4_err);
	g_amp_4kHz->SetPoint(3,amp4,fft4_4);
	g_amp_4kHz->SetPointError(3,amperr4,fft4_4_err);

	new TCanvas();
	g_amp_4kHz->GetXaxis()->SetTitle("Exp amplitude [mV]");
	g_amp_4kHz->GetYaxis()->SetTitle("FFT 4kHz amplitude [arb. u.]");
	g_amp_4kHz->SetMarkerStyle(20);
	g_amp_4kHz->Draw("APL");



	TCanvas* can_g = new TCanvas("can_g","",1800,600);
	can_g->Divide(3,1);
	can_g->cd(1);
	g_kick_blum->Draw("APL");
	can_g->cd(2);
	g_amp_17kHz->Draw("APL");
	can_g->cd(3);
	g_amp_4kHz->Draw("APL");

	gStyle->SetOptStat(0);

	TCanvas* can = new TCanvas("can","",1800,900);
	can->Divide(2,1);
	can->cd(1);
	f1_kick1_exp->SetLineColor(1);
	f2_kick1_exp->SetLineColor(2);
	f3_kick1_exp->SetLineColor(3);
	f4_kick1_exp->SetLineColor(4);
	f1_kick1_exp->Draw("HIST");
	f2_kick1_exp->Draw("HIST SAME");
	f3_kick1_exp->Draw("HIST SAME");
	f4_kick1_exp->Draw("HIST SAME");
	f1_kick1_exp->GetXaxis()->SetRangeUser(-1,1);
	f1_kick1_exp->GetYaxis()->SetRangeUser(-30,60);
	f2_kick1_exp->GetXaxis()->SetRangeUser(-1,1);
	f2_kick1_exp->GetYaxis()->SetRangeUser(-30,60);
	f3_kick1_exp->GetXaxis()->SetRangeUser(-1,1);
	f3_kick1_exp->GetYaxis()->SetRangeUser(-30,60);
	f4_kick1_exp->GetXaxis()->SetRangeUser(-1,1);
	f4_kick1_exp->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	TLegend* leg1 = new TLegend();
	leg1->AddEntry(f1_kick1_exp,Form("Kick %.2f kV",strengths[0]),"L");
	leg1->AddEntry(f2_kick1_exp,Form("Kick %.2f kV",strengths[1]),"L");
	leg1->AddEntry(f3_kick1_exp,Form("Kick %.2f kV",strengths[2]),"L");
	leg1->AddEntry(f4_kick1_exp,Form("Kick %.2f kV",strengths[3]),"L");
	leg1->Draw();
	can->cd(2);
	f1_kick8_exp->SetLineColor(1);
	f2_kick8_exp->SetLineColor(2);
	f3_kick8_exp->SetLineColor(3);
	f4_kick8_exp->SetLineColor(4);
	f1_kick8_exp->Draw("HIST");
	f2_kick8_exp->Draw("HIST SAME");
	f3_kick8_exp->Draw("HIST SAME");
	f4_kick8_exp->Draw("HIST SAME");
	f1_kick8_exp->GetXaxis()->SetRangeUser(-1,1);
	f1_kick8_exp->GetYaxis()->SetRangeUser(-30,60);
	f2_kick8_exp->GetXaxis()->SetRangeUser(-1,1);
	f2_kick8_exp->GetYaxis()->SetRangeUser(-30,60);
	f3_kick8_exp->GetXaxis()->SetRangeUser(-1,1);
	f3_kick8_exp->GetYaxis()->SetRangeUser(-30,60);
	f4_kick8_exp->GetXaxis()->SetRangeUser(-1,1);
	f4_kick8_exp->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();



	new TCanvas();
	f1_kick1_ra->SetLineColor(1);
	f2_kick1_ra->SetLineColor(2);
	f3_kick1_ra->SetLineColor(3);
	f4_kick1_ra->SetLineColor(4);
	f1_kick1_ra->Draw("HIST");
	f2_kick1_ra->Draw("HIST SAME");
	f3_kick1_ra->Draw("HIST SAME");
	f4_kick1_ra->Draw("HIST SAME");
	f1_kick1_ra->GetXaxis()->SetRangeUser(-1,2);
	f1_kick1_ra->GetYaxis()->SetRangeUser(-30,60);
	f2_kick1_ra->GetXaxis()->SetRangeUser(-1,2);
	f2_kick1_ra->GetYaxis()->SetRangeUser(-30,60);
	f3_kick1_ra->GetXaxis()->SetRangeUser(-1,2);
	f3_kick1_ra->GetYaxis()->SetRangeUser(-30,60);
	f4_kick1_ra->GetXaxis()->SetRangeUser(-1,2);
	f4_kick1_ra->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();


	TCanvas* can2 = new TCanvas("can2","",1800,900);
	can2->Divide(2,1);
	can2->cd(1);
	f1_kick1->SetLineColor(1);
	f2_kick1->SetLineColor(2);
	f3_kick1->SetLineColor(3);
	f4_kick1->SetLineColor(4);
	f1_kick1->Draw("HIST");
	f2_kick1->Draw("HIST SAME");
	f3_kick1->Draw("HIST SAME");
	f4_kick1->Draw("HIST SAME");
	f1_kick1->GetXaxis()->SetRangeUser(-1,1);
	f1_kick1->GetYaxis()->SetRangeUser(-30,60);
	f2_kick1->GetXaxis()->SetRangeUser(-1,1);
	f2_kick1->GetYaxis()->SetRangeUser(-30,60);
	f3_kick1->GetXaxis()->SetRangeUser(-1,1);
	f3_kick1->GetYaxis()->SetRangeUser(-30,60);
	f4_kick1->GetXaxis()->SetRangeUser(-1,1);
	f4_kick1->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();
	can2->cd(2);
	f1_kick8->SetLineColor(1);
	f2_kick8->SetLineColor(2);
	f3_kick8->SetLineColor(3);
	f4_kick8->SetLineColor(4);
	f1_kick8->Draw("HIST");
	f2_kick8->Draw("HIST SAME");
	f3_kick8->Draw("HIST SAME");
	f4_kick8->Draw("HIST SAME");
	f1_kick8->GetXaxis()->SetRangeUser(-1,1);
	f1_kick8->GetYaxis()->SetRangeUser(-30,60);
	f2_kick8->GetXaxis()->SetRangeUser(-1,1);
	f2_kick8->GetYaxis()->SetRangeUser(-30,60);
	f3_kick8->GetXaxis()->SetRangeUser(-1,1);
	f3_kick8->GetYaxis()->SetRangeUser(-30,60);
	f4_kick8->GetXaxis()->SetRangeUser(-1,1);
	f4_kick8->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();


	TCanvas* can_FFT = new TCanvas("can_FFT","",1800,900);
	f1_kick1_FFT->SetLineColor(1);
	f2_kick1_FFT->SetLineColor(2);
	f3_kick1_FFT->SetLineColor(3);
	f4_kick1_FFT->SetLineColor(4);
	f1_kick1_FFT->Draw("HIST");
    f2_kick1_FFT->Draw("HIST SAME");
    f3_kick1_FFT->Draw("HIST SAME");
    f4_kick1_FFT->Draw("HIST SAME");
    gPad->SetLogx();
	gPad->BuildLegend();



}
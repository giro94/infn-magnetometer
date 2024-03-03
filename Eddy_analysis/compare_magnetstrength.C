#include "../analysis_tools.C"



void compare_magnetstrength(){

	/*
	//TFile* f1 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct12_H20_nofilter_B0.root");
	TFile* f1 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct13_H20_nofilter_B0.root");
	//TFile* f1 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct14_H25_nofilter_B0.root");
	TFile* f2 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct11_H20_nofilter_B3884.root");
	//TFile* f3 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct12_H20_nofilter_B4595.root");
	TFile* f3 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct13_H20_nofilter_B4595.root");
	TFile* f4 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct10_H20_nofilter_Bfield.root");
	vector<double> strengths = {0, 3884, 4595, 5175};
	*/
	
	//October comparison
	/*
	TFile* f1 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct23_H25_B0_k777.root");
	TFile* f2 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct20_H25_B50_k777.root");
	TFile* f3 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct21_H25_B88_k777.root");
	TFile* f4 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct21_H25_B100_k777.root");
	//TFile* f4 = TFile::Open("fitted_noSNR_kick_analysis_output_SD_R0_eddy_oct22_H25_B100_k777.root");
	vector<double> strengths = {0, 2468, 4595, 5175};
	*/

	//January comparison
	TFile* f1 = TFile::Open("fits/fitted_noSNR_kick_analysis_output_EC_jan21_B3043_H24.root");
	TFile* f2 = TFile::Open("fits/fitted_noSNR_kick_analysis_output_EC_jan20_B3619_H16.root");
	TFile* f3 = TFile::Open("fits/fitted_noSNR_kick_analysis_output_EC_jan19_B4353_H23.root");
	TFile* f4 = TFile::Open("fits/fitted_noSNR_kick_analysis_output_EC_jan19_B5175_H15.root");
	vector<double> strengths = {3043, 3619, 4353, 5175};


	TH1D* f1_kick1 = (TH1D*)f1->Get("trace_kick1_px");
	TH1D* f2_kick1 = (TH1D*)f2->Get("trace_kick1_px");
	TH1D* f3_kick1 = (TH1D*)f3->Get("trace_kick1_px");
	TH1D* f4_kick1 = (TH1D*)f4->Get("trace_kick1_px");

	TH1D* f1_kick8 = (TH1D*)f1->Get("trace_kick8_px");
	TH1D* f2_kick8 = (TH1D*)f2->Get("trace_kick8_px");
	TH1D* f3_kick8 = (TH1D*)f3->Get("trace_kick8_px");
	TH1D* f4_kick8 = (TH1D*)f4->Get("trace_kick8_px");

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

	TH1D* f1_kick1_FFT = (TH1D*)f1->Get("trace_kick1_px_FFT");
	TH1D* f2_kick1_FFT = (TH1D*)f2->Get("trace_kick1_px_FFT");
	TH1D* f3_kick1_FFT = (TH1D*)f3->Get("trace_kick1_px_FFT");
	TH1D* f4_kick1_FFT = (TH1D*)f4->Get("trace_kick1_px_FFT");

	f1_kick1->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[0]));
	f2_kick1->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[1]));
	f3_kick1->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[2]));
	f4_kick1->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[3]));

	f1_kick1_exp->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[0]));
	f2_kick1_exp->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[1]));
	f3_kick1_exp->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[2]));
	f4_kick1_exp->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[3]));

	f1_kick8->SetTitle(Form("Kick 8 (Magnet %.f A)",strengths[0]));
	f2_kick8->SetTitle(Form("Kick 8 (Magnet %.f A)",strengths[1]));
	f3_kick8->SetTitle(Form("Kick 8 (Magnet %.f A)",strengths[2]));
	f4_kick8->SetTitle(Form("Kick 8 (Magnet %.f A)",strengths[3]));

	f1_kick1_FFT->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[0]));
	f2_kick1_FFT->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[1]));
	f3_kick1_FFT->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[2]));
	f4_kick1_FFT->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[3]));


	TGraphErrors* g_blum_amp = new TGraphErrors();
	double amp1 = f1_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr1 = f1_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base1 = f1_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr1 = f1_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak1 = f1_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr1 = f1_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum1 = peak1 - base1;
	double blumerr1 = sqrt(baseerr1*baseerr1 + peakerr1*peakerr1);
	g_blum_amp->SetPoint(g_blum_amp->GetN(),blum1,amp1);
	g_blum_amp->SetPointError(g_blum_amp->GetN()-1,blumerr1,amperr1);
	double amp2 = f2_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr2 = f2_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base2 = f2_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr2 = f2_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak2 = f2_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr2 = f2_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum2 = peak2 - base2;
	double blumerr2 = sqrt(baseerr2*baseerr2 + peakerr2*peakerr2);
	g_blum_amp->SetPoint(g_blum_amp->GetN(),blum2,amp2);
	g_blum_amp->SetPointError(g_blum_amp->GetN()-1,blumerr2,amperr2);
	double amp3 = f3_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr3 = f3_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base3 = f3_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr3 = f3_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak3 = f3_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr3 = f3_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum3 = peak3 - base3;
	double blumerr3 = sqrt(baseerr3*baseerr3 + peakerr3*peakerr3);
	g_blum_amp->SetPoint(g_blum_amp->GetN(),blum3,amp3);
	g_blum_amp->SetPointError(g_blum_amp->GetN()-1,blumerr3,amperr3);
	double amp4 = f4_kick1_exp->GetFunction("f_exp")->GetParameter(0);
	double amperr4 = f4_kick1_exp->GetFunction("f_exp")->GetParError(0);
	double base4 = f4_kick1_exp_res->GetFunction("pol0")->GetParameter(0);
	double baseerr4 = f4_kick1_exp_res->GetFunction("pol0")->GetParError(0);
	double peak4 = f4_kick1_exp_res->GetFunction("f_blumlein")->GetParameter(0);
	double peakerr4 = f4_kick1_exp_res->GetFunction("f_blumlein")->GetParError(0);
	double blum4 = peak4 - base4;
	double blumerr4 = sqrt(baseerr4*baseerr4 + peakerr4*peakerr4);
	g_blum_amp->SetPoint(g_blum_amp->GetN(),blum4,amp4);
	g_blum_amp->SetPointError(g_blum_amp->GetN()-1,blumerr4,amperr4);


	vector<double> norm_factor;
	norm_factor.push_back(50./peak1);
	norm_factor.push_back(50./peak2);
	norm_factor.push_back(50./peak3);
	norm_factor.push_back(50./peak4);


	TH1D* f1_kick1_norm = (TH1D*)f1_kick1_exp->Clone("trace_kick1_norm");
	TH1D* f2_kick1_norm = (TH1D*)f2_kick1_exp->Clone("trace_kick1_norm");
	TH1D* f3_kick1_norm = (TH1D*)f3_kick1_exp->Clone("trace_kick1_norm");
	TH1D* f4_kick1_norm = (TH1D*)f4_kick1_exp->Clone("trace_kick1_norm");
	TH1D* f1_kick8_norm = (TH1D*)f1_kick8_exp->Clone("trace_kick8_norm");
	TH1D* f2_kick8_norm = (TH1D*)f2_kick8_exp->Clone("trace_kick8_norm");
	TH1D* f3_kick8_norm = (TH1D*)f3_kick8_exp->Clone("trace_kick8_norm");
	TH1D* f4_kick8_norm = (TH1D*)f4_kick8_exp->Clone("trace_kick8_norm");

	f1_kick1_norm->Scale(norm_factor[0]);
	f2_kick1_norm->Scale(norm_factor[1]);
	f3_kick1_norm->Scale(norm_factor[2]);
	f4_kick1_norm->Scale(norm_factor[3]);

	f1_kick8_norm->Scale(norm_factor[0]);
	f2_kick8_norm->Scale(norm_factor[1]);
	f3_kick8_norm->Scale(norm_factor[2]);
	f4_kick8_norm->Scale(norm_factor[3]);


	f1_kick1_norm->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[0]));
	f2_kick1_norm->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[1]));
	f3_kick1_norm->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[2]));
	f4_kick1_norm->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[3]));
	f1_kick1_norm->SetLineWidth(2);
	f2_kick1_norm->SetLineWidth(2);
	f3_kick1_norm->SetLineWidth(2);
	f4_kick1_norm->SetLineWidth(2);
	f1_kick1_norm->SetLineColor(1);
	f2_kick1_norm->SetLineColor(2);
	f3_kick1_norm->SetLineColor(3);
	f4_kick1_norm->SetLineColor(4);

	TGraphErrors* g_magnet_amp = new TGraphErrors();
	g_magnet_amp->SetPoint(g_magnet_amp->GetN(),strengths[0],amp1*norm_factor[0]);
	g_magnet_amp->SetPointError(g_magnet_amp->GetN()-1,0,amperr1*norm_factor[0]);
	g_magnet_amp->SetPoint(g_magnet_amp->GetN(),strengths[1],amp2*norm_factor[1]);
	g_magnet_amp->SetPointError(g_magnet_amp->GetN()-1,0,amperr2*norm_factor[1]);
	g_magnet_amp->SetPoint(g_magnet_amp->GetN(),strengths[2],amp3*norm_factor[2]);
	g_magnet_amp->SetPointError(g_magnet_amp->GetN()-1,0,amperr3*norm_factor[2]);
	g_magnet_amp->SetPoint(g_magnet_amp->GetN(),strengths[3],amp4*norm_factor[3]);
	g_magnet_amp->SetPointError(g_magnet_amp->GetN()-1,0,amperr4*norm_factor[3]);
	new TCanvas();
	g_magnet_amp->GetXaxis()->SetTitle("Magnet strength [A]");
	g_magnet_amp->GetYaxis()->SetTitle("Exp amplitude [mV]");
	g_magnet_amp->SetMarkerStyle(20);
	g_magnet_amp->Draw("APL");



	new TCanvas();
	g_blum_amp->GetXaxis()->SetTitle("Blumlein [mV]");
	g_blum_amp->GetYaxis()->SetTitle("Exp amplitude [mV]");
	g_blum_amp->SetMarkerStyle(20);
	g_blum_amp->Draw("APL");

	int bin17 = f1_kick1_FFT->FindBin(17);
	int bin18 = f1_kick1_FFT->FindBin(18);
	TGraphErrors* g_magnet_17kHz = new TGraphErrors();
	double fft1_17_err;
	double fft1_17 = f1_kick1_FFT->IntegralAndError(bin17,bin18,fft1_17_err);
	g_magnet_17kHz->SetPoint(g_magnet_17kHz->GetN(),strengths[0],fft1_17);
	g_magnet_17kHz->SetPointError(g_magnet_17kHz->GetN()-1,0,fft1_17_err);
	double fft2_17_err;
	double fft2_17 = f2_kick1_FFT->IntegralAndError(bin17,bin18,fft2_17_err);
	g_magnet_17kHz->SetPoint(g_magnet_17kHz->GetN(),strengths[1],fft2_17);
	g_magnet_17kHz->SetPointError(g_magnet_17kHz->GetN()-1,0,fft2_17_err);
	double fft3_17_err;
	double fft3_17 = f3_kick1_FFT->IntegralAndError(bin17,bin18,fft3_17_err);
	g_magnet_17kHz->SetPoint(g_magnet_17kHz->GetN(),strengths[2],fft3_17);
	g_magnet_17kHz->SetPointError(g_magnet_17kHz->GetN()-1,0,fft3_17_err);
	double fft4_17_err;
	double fft4_17 = f4_kick1_FFT->IntegralAndError(bin17,bin18,fft4_17_err);
	g_magnet_17kHz->SetPoint(g_magnet_17kHz->GetN(),strengths[3],fft4_17);
	g_magnet_17kHz->SetPointError(g_magnet_17kHz->GetN()-1,0,fft4_17_err);

	new TCanvas();
	g_magnet_17kHz->GetXaxis()->SetTitle("Magnet strength [A]");
	g_magnet_17kHz->GetYaxis()->SetTitle("FFT 17kHz amplitude [arb. u.]");
	g_magnet_17kHz->SetMarkerStyle(20);
	g_magnet_17kHz->Draw("APL");


	int bin4 = f1_kick1_FFT->FindBin(4);
	int bin4p5 = f1_kick1_FFT->FindBin(4.5);
	TGraphErrors* g_magnet_4kHz = new TGraphErrors();
	double fft1_4_err;
	double fft1_4 = f1_kick1_FFT->IntegralAndError(bin4,bin18,fft1_4_err);
	g_magnet_4kHz->SetPoint(g_magnet_4kHz->GetN(),strengths[0],fft1_4);
	g_magnet_4kHz->SetPointError(g_magnet_4kHz->GetN()-1,0,fft1_4_err);
	double fft2_4_err;
	double fft2_4 = f2_kick1_FFT->IntegralAndError(bin4,bin18,fft2_4_err);
	g_magnet_4kHz->SetPoint(g_magnet_4kHz->GetN(),strengths[1],fft2_4);
	g_magnet_4kHz->SetPointError(g_magnet_4kHz->GetN()-1,0,fft2_4_err);
	double fft3_4_err;
	double fft3_4 = f3_kick1_FFT->IntegralAndError(bin4,bin18,fft3_4_err);
	g_magnet_4kHz->SetPoint(g_magnet_4kHz->GetN(),strengths[2],fft3_4);
	g_magnet_4kHz->SetPointError(g_magnet_4kHz->GetN()-1,0,fft3_4_err);
	double fft4_4_err;
	double fft4_4 = f4_kick1_FFT->IntegralAndError(bin4,bin18,fft4_4_err);
	g_magnet_4kHz->SetPoint(g_magnet_4kHz->GetN(),strengths[3],fft4_4);
	g_magnet_4kHz->SetPointError(g_magnet_4kHz->GetN()-1,0,fft4_4_err);

	new TCanvas();
	g_magnet_4kHz->GetXaxis()->SetTitle("Magnet strength [A]");
	g_magnet_4kHz->GetYaxis()->SetTitle("FFT 4kHz amplitude [arb. u.]");
	g_magnet_4kHz->SetMarkerStyle(20);
	g_magnet_4kHz->Draw("APL");

	TCanvas* can_g = new TCanvas("can_g","",1800,600);
	can_g->Divide(3,1);
	can_g->cd(1);
	g_magnet_amp->Draw("APL");
	can_g->cd(2);
	g_magnet_17kHz->Draw("APL");
	can_g->cd(3);
	g_magnet_4kHz->Draw("APL");

	gStyle->SetOptStat(0);

	TCanvas* can = new TCanvas("can","",1800,900);
	can->Divide(2,1);
	can->cd(1);
	f1_kick1_exp->SetLineWidth(2);
	f2_kick1_exp->SetLineWidth(2);
	f3_kick1_exp->SetLineWidth(2);
	f4_kick1_exp->SetLineWidth(2);
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
	leg1->AddEntry(f1_kick1_exp,Form("Magnet %.1f A",strengths[0]),"L");
	leg1->AddEntry(f2_kick1_exp,Form("Magnet %.1f A",strengths[1]),"L");
	leg1->AddEntry(f3_kick1_exp,Form("Magnet %.1f A",strengths[2]),"L");
	leg1->AddEntry(f4_kick1_exp,Form("Magnet %.1f A",strengths[3]),"L");
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

	TCanvas* can3 = new TCanvas("can3","",1800,900);
	can3->Divide(2,1);
	can3->cd(1);
	f1_kick1_norm->SetLineColor(1);
	f2_kick1_norm->SetLineColor(2);
	f3_kick1_norm->SetLineColor(3);
	f4_kick1_norm->SetLineColor(4);
	f1_kick1_norm->Draw("HIST");
	f2_kick1_norm->Draw("HIST SAME");
	f3_kick1_norm->Draw("HIST SAME");
	f4_kick1_norm->Draw("HIST SAME");
	f1_kick1_norm->GetXaxis()->SetRangeUser(-1,1);
	f1_kick1_norm->GetYaxis()->SetRangeUser(-30,60);
	f2_kick1_norm->GetXaxis()->SetRangeUser(-1,1);
	f2_kick1_norm->GetYaxis()->SetRangeUser(-30,60);
	f3_kick1_norm->GetXaxis()->SetRangeUser(-1,1);
	f3_kick1_norm->GetYaxis()->SetRangeUser(-30,60);
	f4_kick1_norm->GetXaxis()->SetRangeUser(-1,1);
	f4_kick1_norm->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();
	can3->cd(2);
	f1_kick8_norm->SetLineColor(1);
	f2_kick8_norm->SetLineColor(2);
	f3_kick8_norm->SetLineColor(3);
	f4_kick8_norm->SetLineColor(4);
	f1_kick8_norm->Draw("HIST");
	f2_kick8_norm->Draw("HIST SAME");
	f3_kick8_norm->Draw("HIST SAME");
	f4_kick8_norm->Draw("HIST SAME");
	f1_kick8_norm->GetXaxis()->SetRangeUser(-1,1);
	f1_kick8_norm->GetYaxis()->SetRangeUser(-30,60);
	f2_kick8_norm->GetXaxis()->SetRangeUser(-1,1);
	f2_kick8_norm->GetYaxis()->SetRangeUser(-30,60);
	f3_kick8_norm->GetXaxis()->SetRangeUser(-1,1);
	f3_kick8_norm->GetYaxis()->SetRangeUser(-30,60);
	f4_kick8_norm->GetXaxis()->SetRangeUser(-1,1);
	f4_kick8_norm->GetYaxis()->SetRangeUser(-30,60);
	gPad->SetGridy();
	gPad->BuildLegend();



}
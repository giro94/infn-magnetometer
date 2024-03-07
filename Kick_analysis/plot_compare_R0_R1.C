#include "../analysis_tools.C"

void plot_compare_R0_R1(){

	TFile* f0 = TFile::Open("analysis_FD_R0_oct9_H0_Bfield.root");
	TFile* f1 = TFile::Open("analysis_FD_R1_oct8_H0_Bfield.root");

	TFile* f0_blum = TFile::Open("../Eddy_analysis/analysis/analysis_FD_R0_blum_oct9_H22p5_Bfield.root");
	TFile* f1_blum = TFile::Open("../Eddy_analysis/analysis/analysis_FD_R1_blum_oct8_H0_Bfield.root");

	TFile* f0_ramp = TFile::Open("../Ramp_up_analysis/rampup_analysis_output_FD_R0_ramp_oct9_H22p5.root");
	TFile* f1_ramp = TFile::Open("../Ramp_up_analysis/rampup_analysis_output_FD_R1_ramp_oct8_H0.root");

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;



	TGraph* h1_kick1_R0 = (TGraph*)f0->Get("kick1_rep1");
	TGraph* h1_lastkick_R0 = (TGraph*)f0->Get("kick8_rep4");
	TGraph*	h1_kick1_R0_norm = (TGraph*)f0->Get("normalized_kick_1");
	TGraph* h1_kick1_R1 = (TGraph*)f1->Get("kick1_rep1");
	TGraph* h1_lastkick_R1 = (TGraph*)f1->Get("kick8_rep5");
	TGraph*	h1_kick1_R1_norm = (TGraph*)f1->Get("normalized_kick_1");
	h1_kick1_R0_norm->SetTitle("Kick R0 (normalized)");
	h1_kick1_R1_norm->SetTitle("Kick R1 (normalized)");



	TH1D* h1_kick1_R0_blum = ((TProfile*)f0_blum->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R1_blum = ((TProfile*)f1_blum->Get("trace_kick1"))->ProjectionX();
	h1_kick1_R0_blum->Scale(-1);

	TH1D* h1_kick1_R0_blum_ra = smoothing(h1_kick1_R0_blum,"");
	TH1D* h1_kick1_R1_blum_ra = smoothing(h1_kick1_R1_blum,"");
	h1_kick1_R0_blum_ra->SetTitle("Blumlein R0");
	h1_kick1_R1_blum_ra->SetTitle("Blumlein R1");








	new TCanvas();
	h1_kick1_R0_norm->SetLineWidth(2);
	h1_kick1_R1_norm->SetLineWidth(2);
	h1_kick1_R0_norm->SetLineColor(kBlack);
	h1_kick1_R1_norm->SetLineColor(kRed);
	h1_kick1_R0_norm->Draw("AL");
	h1_kick1_R1_norm->Draw("L");
	gPad->BuildLegend();




	new TCanvas();
	h1_kick1_R0_blum_ra->SetLineWidth(2);
	h1_kick1_R1_blum_ra->SetLineWidth(2);
	h1_kick1_R0_blum_ra->SetLineColor(kBlack);
	h1_kick1_R1_blum_ra->SetLineColor(kRed);
	h1_kick1_R0_blum_ra->Draw("HIST");
	h1_kick1_R1_blum_ra->Draw("HIST SAME");
	gPad->BuildLegend();

	TF1* f_baseline = new TF1("f_baseline","[0]",-1,-0.55);
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.4,-0.2);
	f_blumlein->SetParameters(0.5,0.3,-1000.0);

	TFitResultPtr fit_R0_blum = h1_kick1_R0_blum_ra->Fit("f_blumlein","QSR","");
	TFitResultPtr fit_R1_blum = h1_kick1_R1_blum_ra->Fit("f_blumlein","QSR","");

	TFitResultPtr fit_R0_base = h1_kick1_R0_blum_ra->Fit("f_baseline","QSR","");
	TFitResultPtr fit_R1_base = h1_kick1_R1_blum_ra->Fit("f_baseline","QSR","");

	double R0_blum = fit_R0_blum->Parameter(0) - fit_R0_base->Parameter(0);
	double R1_blum = fit_R1_blum->Parameter(0) - fit_R1_base->Parameter(0);
	cout<<"Blum R0 : "<<R0_blum<<"\n";
	cout<<"Blum R1 : "<<R1_blum<<"\n";

	for (int i=0; i<h1_lastkick_R0->GetN(); i++){
		h1_lastkick_R0->SetPoint(i,h1_lastkick_R0->GetPointX(i),h1_lastkick_R0->GetPointY(i)*129/R0_blum);
	}
	for (int i=0; i<h1_lastkick_R1->GetN(); i++){
		h1_lastkick_R1->SetPoint(i,h1_lastkick_R1->GetPointX(i),h1_lastkick_R1->GetPointY(i)*129/R1_blum);
	}

	new TCanvas();
	new TCanvas();
	h1_lastkick_R0->SetLineWidth(2);
	h1_lastkick_R1->SetLineWidth(2);
	h1_lastkick_R0->SetLineColor(kBlack);
	h1_lastkick_R1->SetLineColor(kRed);
	h1_lastkick_R0->Draw("AL");
	h1_lastkick_R1->Draw("L");
	gPad->BuildLegend();






}
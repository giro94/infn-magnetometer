#include "../analysis_tools.C"


void compare_quadrant(){


	TFile* f1 = TFile::Open("output_Quadrant_R1_Bfield_oct7_H0_k000.root");
	TFile* f2 = TFile::Open("output_Quadrant_R1_Bfield_oct7_H0_k194.root");
	TFile* f3 = TFile::Open("output_Quadrant_R1_Bfield_oct7_H0_k389.root");
	TFile* f4 = TFile::Open("output_Quadrant_R1_Bfield_oct7_H0_k583.root");
	TFile* f5 = TFile::Open("output_Quadrant_R1_Bfield_oct7_H0_k777.root");
	vector<double> strengths = {0.00, 1.94, 3.89, 5.83, 7.77};

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1x = ((TProfile*)f1->Get("p_traceX"))->ProjectionX("h1x");
	TH1D* h2x = ((TProfile*)f2->Get("p_traceX"))->ProjectionX("h2x");
	TH1D* h3x = ((TProfile*)f3->Get("p_traceX"))->ProjectionX("h3x");
	TH1D* h4x = ((TProfile*)f4->Get("p_traceX"))->ProjectionX("h4x");
	TH1D* h5x = ((TProfile*)f5->Get("p_traceX"))->ProjectionX("h5x");
	TH1D* h1y = ((TProfile*)f1->Get("p_traceY"))->ProjectionX("h1y");
	TH1D* h2y = ((TProfile*)f2->Get("p_traceY"))->ProjectionX("h2y");
	TH1D* h3y = ((TProfile*)f3->Get("p_traceY"))->ProjectionX("h3y");
	TH1D* h4y = ((TProfile*)f4->Get("p_traceY"))->ProjectionX("h4y");
	TH1D* h5y = ((TProfile*)f5->Get("p_traceY"))->ProjectionX("h5y");

	TH1D* h1x_ra = smoothing(h1x,"h1x_ra");
	TH1D* h2x_ra = smoothing(h2x,"h2x_ra");
	TH1D* h3x_ra = smoothing(h3x,"h3x_ra");
	TH1D* h4x_ra = smoothing(h4x,"h4x_ra");
	TH1D* h5x_ra = smoothing(h5x,"h5x_ra");
	TH1D* h1y_ra = smoothing(h1y,"h1y_ra");
	TH1D* h2y_ra = smoothing(h2y,"h2y_ra");
	TH1D* h3y_ra = smoothing(h3y,"h3y_ra");
	TH1D* h4y_ra = smoothing(h4y,"h4y_ra");
	TH1D* h5y_ra = smoothing(h5y,"h5y_ra");


	vector<TLine*> lines;
	for (int i=0; i<8; i++){
		lines.push_back(new TLine(5+i*10,-10,5+i*10,10));
		lines.back()->SetLineStyle(kDashed);
		lines.back()->SetLineColor(kBlack);
		lines.back()->SetLineWidth(2);
	}

	gStyle->SetOptStat(0);
	new TCanvas("","",1200,600);
	h1x_ra->GetYaxis()->SetRangeUser(-10,10);
	h1x_ra->SetLineWidth(2);
	h2x_ra->SetLineWidth(2);
	h3x_ra->SetLineWidth(2);
	h4x_ra->SetLineWidth(2);
	h5x_ra->SetLineWidth(2);
	h1x_ra->SetLineColor(kOrange);
	h2x_ra->SetLineColor(kBlack);
	h3x_ra->SetLineColor(kRed);
	h4x_ra->SetLineColor(kGreen);
	h5x_ra->SetLineColor(kBlue);
	h1x_ra->GetYaxis()->SetTitle("X displacement [arb.u.]");
	h2x_ra->GetYaxis()->SetTitle("X displacement [arb.u.]");
	h3x_ra->GetYaxis()->SetTitle("X displacement [arb.u.]");
	h4x_ra->GetYaxis()->SetTitle("X displacement [arb.u.]");
	h5x_ra->GetYaxis()->SetTitle("X displacement [arb.u.]");
	h1x_ra->SetTitle(Form("Kick %.f %%",100.*strengths[0]/7.77));
	h2x_ra->SetTitle(Form("Kick %.f %%",100.*strengths[1]/7.77));
	h3x_ra->SetTitle(Form("Kick %.f %%",100.*strengths[2]/7.77));
	h4x_ra->SetTitle(Form("Kick %.f %%",100.*strengths[3]/7.77));
	h5x_ra->SetTitle(Form("Kick %.f %%",100.*strengths[4]/7.77));
	h1x_ra->Draw("HIST");
	h2x_ra->Draw("HIST SAME");
	h3x_ra->Draw("HIST SAME");
	h4x_ra->Draw("HIST SAME");
	h5x_ra->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend();
	for (auto line : lines) line->Draw("SAME");

	new TCanvas("","",1200,600);
	h1y_ra->GetYaxis()->SetRangeUser(-10,10);
	h1y_ra->SetLineWidth(2);
	h2y_ra->SetLineWidth(2);
	h3y_ra->SetLineWidth(2);
	h4y_ra->SetLineWidth(2);
	h5y_ra->SetLineWidth(2);
	h1y_ra->SetLineColor(kOrange);
	h2y_ra->SetLineColor(kBlack);
	h3y_ra->SetLineColor(kRed);
	h4y_ra->SetLineColor(kGreen);
	h5y_ra->SetLineColor(kBlue);
	h1y_ra->GetYaxis()->SetTitle("Y displacement [arb.u.]");
	h2y_ra->GetYaxis()->SetTitle("Y displacement [arb.u.]");
	h3y_ra->GetYaxis()->SetTitle("Y displacement [arb.u.]");
	h4y_ra->GetYaxis()->SetTitle("Y displacement [arb.u.]");
	h5y_ra->GetYaxis()->SetTitle("Y displacement [arb.u.]");
	h1y_ra->SetTitle(Form("Kick %.f %%",100.*strengths[0]/7.77));
	h2y_ra->SetTitle(Form("Kick %.f %%",100.*strengths[1]/7.77));
	h3y_ra->SetTitle(Form("Kick %.f %%",100.*strengths[2]/7.77));
	h4y_ra->SetTitle(Form("Kick %.f %%",100.*strengths[3]/7.77));
	h5y_ra->SetTitle(Form("Kick %.f %%",100.*strengths[4]/7.77));
	h1y_ra->Draw("HIST");
	h2y_ra->Draw("HIST SAME");
	h3y_ra->Draw("HIST SAME");
	h4y_ra->Draw("HIST SAME");
	h5y_ra->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend();
	for (auto line : lines) line->Draw("SAME");


}
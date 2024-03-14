#include "../analysis_tools.C"

void plot_compare_R0_R1(){

	TFile* f0 = TFile::Open("output_FD_R0_oct9_H0_Bfield.root");
	TFile* f1 = TFile::Open("output_FD_R1_oct8_H0_Bfield.root");

	TFile* f0_blum = TFile::Open("analysis_FD_R0_blum_oct9_H22p5_Bfield.root");
	TFile* f1_blum = TFile::Open("analysis_FD_R1_blum_oct8_H0_Bfield.root");

	TFile* f0_ramp = TFile::Open("output_FD_R0_ramp_oct9_H22p5.root");
	TFile* f1_ramp = TFile::Open("output_FD_R1_ramp_oct8_H0.root");

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;



	TGraph* h1_kick1_R0 = (TGraph*)f0->Get("kick1_rep1");
	TGraph* h1_lastkick_R0 = (TGraph*)f0->Get("kick8_rep4");
	TGraph*	h1_kick1_R0_norm = (TGraph*)f0->Get("normalized_kick_1");
	TGraph* h1_kick1_R1 = (TGraph*)f1->Get("kick1_rep1");
	TGraph* h1_lastkick_R1 = (TGraph*)f1->Get("kick8_rep5");
	TGraph*	h1_kick1_R1_norm = (TGraph*)f1->Get("normalized_kick_1");
	h1_kick1_R0_norm->SetTitle("Kick R0 (normalized)");
	h1_kick1_R1_norm->SetTitle("Kick R1 (normalized)");

	TGraph* h1_lastkick_R0_calib = new TGraph();
	TGraph* h1_lastkick_R1_calib = new TGraph(); 

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
	h1_kick1_R0_norm->SetLineColor(8);
	h1_kick1_R1_norm->SetLineColor(kBlue);
	h1_kick1_R0_norm->Draw("AL");
	h1_kick1_R1_norm->Draw("L");
	gPad->BuildLegend();




	new TCanvas();
	h1_kick1_R0_blum_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_R0_blum_ra->GetYaxis()->SetRangeUser(-1,1);
	h1_kick1_R1_blum_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_R1_blum_ra->GetYaxis()->SetRangeUser(-1,1);
	h1_kick1_R0_blum_ra->SetLineWidth(2);
	h1_kick1_R1_blum_ra->SetLineWidth(2);
	h1_kick1_R0_blum_ra->SetLineColor(8);
	h1_kick1_R1_blum_ra->SetLineColor(kBlue);
	h1_kick1_R0_blum_ra->Draw("HIST");
	h1_kick1_R1_blum_ra->Draw("HIST SAME");
	gPad->BuildLegend();

	TF1* f_baseline = new TF1("f_baseline","[0]",-1,-0.55);
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.4,-0.2);
	f_blumlein->SetParameters(0.5,-0.3,-1000.0);
	TF1* f_peak = new TF1("f_peak","[0]+[2]*(x-[1])*(x-[1])",-0.05,0.05);
	f_peak->SetParameters(1,0,-1000.0);

	TFitResultPtr fit_R0_blum = h1_kick1_R0_blum_ra->Fit("f_blumlein","QS","",-0.4,-0.2);
	f_blumlein->DrawCopy("SAME");
	TFitResultPtr fit_R1_blum = h1_kick1_R1_blum_ra->Fit("f_blumlein","QS","",-0.4,-0.2);
	f_blumlein->DrawCopy("SAME");

	TFitResultPtr fit_R0_base = h1_kick1_R0_blum_ra->Fit("f_baseline","QS","",-1,-0.55);
	f_baseline->DrawCopy("SAME");
	TFitResultPtr fit_R1_base = h1_kick1_R1_blum_ra->Fit("f_baseline","QS","",-1,-0.55);
	f_baseline->DrawCopy("SAME");

	double R0_blum = fit_R0_blum->Parameter(0) - fit_R0_base->Parameter(0);
	double R1_blum = fit_R1_blum->Parameter(0) - fit_R1_base->Parameter(0);
	double R0_blum_err = sqrt(fit_R0_blum->ParError(0)*fit_R0_blum->ParError(0) + fit_R0_base->ParError(0)*fit_R0_base->ParError(0));
	double R1_blum_err = sqrt(fit_R1_blum->ParError(0)*fit_R1_blum->ParError(0) + fit_R1_base->ParError(0)*fit_R1_base->ParError(0));
	cout<<"Blum R0 : "<<R0_blum<<" +- "<<R0_blum_err<<"\n";
	cout<<"Blum R1 : "<<R1_blum<<" +- "<<R1_blum_err<<"\n";



	new TCanvas();
	h1_lastkick_R0->GetXaxis()->SetRangeUser(-1,1);
	h1_lastkick_R1->GetXaxis()->SetRangeUser(-1,1);
	h1_lastkick_R0->GetYaxis()->SetRangeUser(-0.5,1.2);
	h1_lastkick_R1->GetYaxis()->SetRangeUser(-0.5,1.2);
	h1_lastkick_R0->SetLineWidth(2);
	h1_lastkick_R1->SetLineWidth(2);
	h1_lastkick_R0->SetLineColor(8);
	h1_lastkick_R1->SetLineColor(kBlue);
	h1_lastkick_R0->Draw("AL");
	h1_lastkick_R1->Draw("L");
	gPad->BuildLegend();

	TFitResultPtr fit_R0_peak = h1_lastkick_R0->Fit("f_peak","QS","",-0.05,0.05);
	f_peak->DrawCopy("SAME");
	TFitResultPtr fit_R1_peak = h1_lastkick_R1->Fit("f_peak","QS","",-0.05,0.05);
	f_peak->DrawCopy("SAME");
	TFitResultPtr fit_R0_peakbase = h1_lastkick_R0->Fit("f_baseline","QS","",-1,-0.2);
	f_baseline->DrawCopy("SAME");
	TFitResultPtr fit_R1_peakbase = h1_lastkick_R1->Fit("f_baseline","QS","",-1,-0.2);
	f_baseline->DrawCopy("SAME");
	double R0_peak = fit_R0_peak->Parameter(0) - fit_R0_peakbase->Parameter(0);
	double R1_peak = fit_R1_peak->Parameter(0) - fit_R1_peakbase->Parameter(0);
	double R0_peak_err = sqrt(fit_R0_peak->ParError(0)*fit_R0_peak->ParError(0) + fit_R0_peakbase->ParError(0)*fit_R0_peakbase->ParError(0));
	double R1_peak_err = sqrt(fit_R1_peak->ParError(0)*fit_R1_peak->ParError(0) + fit_R1_peakbase->ParError(0)*fit_R1_peakbase->ParError(0));
	cout<<"Peak R0 : "<<R0_peak<<" +- "<<R0_peak_err<<"\n";
	cout<<"Peak R1 : "<<R1_peak<<" +- "<<R1_peak_err<<"\n";

	for (int i=0; i<h1_lastkick_R0->GetN(); i++){
		h1_lastkick_R0_calib->SetPoint(i,h1_lastkick_R0->GetPointX(i),h1_lastkick_R0->GetPointY(i)*129/R0_blum);
	}
	for (int i=0; i<h1_lastkick_R1->GetN(); i++){
		h1_lastkick_R1_calib->SetPoint(i,h1_lastkick_R1->GetPointX(i),h1_lastkick_R1->GetPointY(i)*129/R1_blum);
	}





	new TCanvas();
	h1_lastkick_R0_calib->SetTitle("Kick R0 calibrated");
	h1_lastkick_R1_calib->SetTitle("Kick R1 calibrated");
	h1_lastkick_R0_calib->GetXaxis()->SetTitle("Time [#mus]");
	h1_lastkick_R1_calib->GetXaxis()->SetTitle("Time [#mus]");
	h1_lastkick_R0_calib->GetYaxis()->SetTitle("B field [G]");
	h1_lastkick_R1_calib->GetYaxis()->SetTitle("B field [G]");
	h1_lastkick_R0_calib->GetXaxis()->SetRangeUser(-0.6,3.5);
	h1_lastkick_R1_calib->GetXaxis()->SetRangeUser(-0.6,3.5);
	h1_lastkick_R0_calib->SetLineWidth(2);
	h1_lastkick_R1_calib->SetLineWidth(2);
	h1_lastkick_R0_calib->SetLineColor(8);
	h1_lastkick_R1_calib->SetLineColor(kBlue);
	h1_lastkick_R0_calib->Draw("AL");
	h1_lastkick_R1_calib->Draw("L");
	gPad->BuildLegend();


	double R0_peak_blum_ratio = R0_peak*1000/R0_blum;
	double R1_peak_blum_ratio = R1_peak*1000/R1_blum;
	double R0_peak_blum_ratio_err = R0_peak_blum_ratio*sqrt(R0_peak_err*R0_peak_err + R0_blum_err*R0_blum_err);
	double R1_peak_blum_ratio_err = R1_peak_blum_ratio*sqrt(R1_peak_err*R1_peak_err + R1_blum_err*R1_blum_err);

	TH1D* g_ratio = new TH1D("g_ratio","",2,0,2);
	g_ratio->SetBinContent(1,R0_peak_blum_ratio);
	g_ratio->SetBinError(1,R0_peak_blum_ratio_err);
	g_ratio->SetBinContent(2,R1_peak_blum_ratio);
	g_ratio->SetBinError(2,R1_peak_blum_ratio_err);

	new TCanvas();
	g_ratio->SetMarkerStyle(20);
	g_ratio->GetXaxis()->SetBinLabel(1,"R0");
	g_ratio->GetXaxis()->SetBinLabel(2,"R1");
	g_ratio->GetYaxis()->SetTitle("Kick/Blumlein [mV/mV]");
	g_ratio->Draw("EX0");
	gPad->SetGridy();
}
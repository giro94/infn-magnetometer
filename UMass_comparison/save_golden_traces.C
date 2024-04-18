#include "../analysis_tools.C"

void save_golden_traces(){

	//R0

	TFile* f1_p = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f1_n = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan28_B5173_H25_Q00.root");
	TFile* f1_0 = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");
	TH1D* h1_kick1_R0_p = ((TProfile*)f1_p->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_n = ((TProfile*)f1_n->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_0 = ((TProfile*)f1_0->Get("trace_kick1"))->ProjectionX();


	double blum_norm_x = -0.323;
	double blum_norm_y_R0 = 129.;
	double blum_norm_y_R1 = 157.;
	double blum_norm_x_trace = 4.769;
	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1_kick1_R0 = (TH1D*)h1_kick1_R0_p->Clone("h1_kick1_R0");
	h1_kick1_R0->Scale(0.5);
	h1_kick1_R0->Add(h1_kick1_R0_n,-0.5);
	h1_kick1_R0->Add(h1_kick1_R0_0,0.2);

	TH1D* h1_kick1_R0_zero = (TH1D*)h1_kick1_R0->Clone("h1_kick1_R0_zero");
	cleanTrace(h1_kick1_R0_zero,-200);
	TH1D* h1_kick1_R0_zero_ra = smoothing(h1_kick1_R0_zero,"");


	TH1D* h1_kick1_R0_ra = smoothing(h1_kick1_R0,"");
	double r0_norm = blum_norm_y_R0/h1_kick1_R0_ra->Interpolate(blum_norm_x);
	h1_kick1_R0_ra->Scale(r0_norm);
	h1_kick1_R0->Scale(r0_norm);
	h1_kick1_R0_zero->Scale(r0_norm);
	h1_kick1_R0_zero_ra->Scale(r0_norm);


	ofstream fout_R0;
	fout_R0.open("INFN_EC_R0_Bon.csv");
	fout_R0<<"Time [ms],Trace [mG]\n";
	for (int bn=1; bn<=h1_kick1_R0_zero_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_zero_ra->GetBinCenter(bn);
		double y = h1_kick1_R0_zero_ra->GetBinContent(bn);
		if (x<-1||x>2) continue;
		fout_R0<<x<<","<<y<<"\n";
	}
	fout_R0.close();

	ofstream fout_R0_raw;
	fout_R0_raw.open("INFN_EC_R0_Bon_raw.csv");
	fout_R0_raw<<"Time [ms],Trace [mG]\n";
	for (int bn=1; bn<=h1_kick1_R0->GetNbinsX(); bn++){
		double x = h1_kick1_R0->GetBinCenter(bn);
		double y = h1_kick1_R0->GetBinContent(bn);
		if (x<-1||x>2) continue;
		fout_R0_raw<<x<<","<<y<<"\n";
	}
	fout_R0_raw.close();



	//R1

	TFile* f_R1 = TFile::Open("../Eddy_analysis/analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root");
	TH1D* h1_kick1_R1 = ((TProfile*)f_R1->Get("trace_kick1"))->ProjectionX();


	TH1D* h1_kick1_R1_zero = (TH1D*)h1_kick1_R1->Clone("h1_kick1_R1_zero");
	cleanTrace(h1_kick1_R1_zero,-200);
	TH1D* h1_kick1_R1_zero_ra = smoothing(h1_kick1_R1_zero,"");

	TH1D* h1_kick1_R1_ra = smoothing(h1_kick1_R1,"");
	double r1_norm = blum_norm_y_R1/h1_kick1_R1_ra->Interpolate(blum_norm_x);
	h1_kick1_R1_ra->Scale(r1_norm);
	h1_kick1_R1_zero_ra->Scale(r1_norm);



	ofstream fout_R1;
	fout_R1.open("INFN_EC_R1_Bon.csv");
	fout_R1<<"Time [ms],Trace [mG]\n";
	for (int bn=1; bn<=h1_kick1_R1_zero_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R1_zero_ra->GetBinCenter(bn);
		double y = h1_kick1_R1_zero_ra->GetBinContent(bn);
		if (x<-1||x>2) continue;
		fout_R1<<x<<","<<y<<"\n";
	}
	fout_R1.close();


	//Fast kick


	TFile* f0 = TFile::Open("../Kick_analysis/output_FD_R0_oct9_H0_Bfield.root");
	TFile* f1 = TFile::Open("../Kick_analysis/output_FD_R1_oct8_H0_Bfield.root");
	TGraph*	h1_kick1_R0_norm = (TGraph*)f0->Get("normalized_kick_1");
	TGraph*	h1_kick1_R1_norm = (TGraph*)f1->Get("normalized_kick_1");


	ofstream fout_kick_R0;
	fout_kick_R0.open("INFN_kick_R0_Bon.csv");
	fout_kick_R0<<"Time [us],Trace [arb. u.]\n";
	for (int bn=0; bn<h1_kick1_R0_norm->GetN(); bn++){
		double x = h1_kick1_R0_norm->GetPointX(bn);
		double y = h1_kick1_R0_norm->GetPointY(bn);
		fout_kick_R0<<x<<","<<y<<"\n";
	}
	fout_kick_R0.close();

	ofstream fout_kick_R1;
	fout_kick_R1.open("INFN_kick_R1_Bon.csv");
	fout_kick_R1<<"Time [us],Trace [arb. u.]\n";
	for (int bn=0; bn<h1_kick1_R1_norm->GetN(); bn++){
		double x = h1_kick1_R1_norm->GetPointX(bn);
		double y = h1_kick1_R1_norm->GetPointY(bn);
		fout_kick_R1<<x<<","<<y<<"\n";
	}
	fout_kick_R1.close();


	new TCanvas();
	h1_kick1_R0->SetLineWidth(2);
	h1_kick1_R0_zero->SetLineWidth(2);
	h1_kick1_R0_ra->SetLineWidth(2);
	h1_kick1_R0_zero_ra->SetLineWidth(2);
	h1_kick1_R0->SetLineColor(kBlack);
	h1_kick1_R0_zero->SetLineColor(kBlack);
	h1_kick1_R0_ra->SetLineColor(kBlue);
	h1_kick1_R0_zero_ra->SetLineColor(kRed);
	h1_kick1_R0->Draw("HIST L");
	h1_kick1_R0_ra->Draw("HIST L SAME");
	h1_kick1_R0_zero->Draw("HIST L SAME");
	h1_kick1_R0_zero_ra->Draw("HIST L SAME");


	new TCanvas();
	h1_kick1_R0_zero_ra->Draw("HIST L");
	h1_kick1_R1_zero_ra->Draw("HIST L SAME");

}
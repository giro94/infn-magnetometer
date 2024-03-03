#include "../analysis_tools.C"

void plot_compare_magnetstrength(){


	TFile* f1 = TFile::Open("analysis/analysis_EC_jan21_B3043_H24.root");
	TFile* f2 = TFile::Open("analysis/analysis_EC_jan20_B3619_H16.root");
	TFile* f3 = TFile::Open("analysis/analysis_EC_jan19_B4353_H23.root");
	TFile* f4 = TFile::Open("analysis/analysis_EC_jan19_B5175_H15.root");
	vector<double> strengths = {3043, 3619, 4353, 5175};


	TH1D* f1_trace = ((TProfile*)f1->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* f2_trace = ((TProfile*)f2->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* f3_trace = ((TProfile*)f3->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* f4_trace = ((TProfile*)f4->Get("trace_kick1_calibrated"))->ProjectionX();


	TH1D* f1_trace_ra = f1_trace; //runningAverage_5_10_15(f1_trace);
	TH1D* f2_trace_ra = f2_trace; //runningAverage_5_10_15(f2_trace);
	TH1D* f3_trace_ra = f3_trace; //runningAverage_5_10_15(f3_trace);
	TH1D* f4_trace_ra = f4_trace; //runningAverage_5_10_15(f4_trace);


	f1_trace_ra->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[0]));
	f2_trace_ra->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[1]));
	f3_trace_ra->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[2]));
	f4_trace_ra->SetTitle(Form("Kick 1 (Magnet %.f A)",strengths[3]));
	f1_trace_ra->GetXaxis()->SetRangeUser(-1,2);
	f2_trace_ra->GetXaxis()->SetRangeUser(-1,2);
	f3_trace_ra->GetXaxis()->SetRangeUser(-1,2);
	f4_trace_ra->GetXaxis()->SetRangeUser(-1,2);
	f1_trace_ra->GetYaxis()->SetRangeUser(-200,200);
	f2_trace_ra->GetYaxis()->SetRangeUser(-200,200);
	f3_trace_ra->GetYaxis()->SetRangeUser(-200,200);
	f4_trace_ra->GetYaxis()->SetRangeUser(-200,200);
	f1_trace_ra->SetLineWidth(2);
	f2_trace_ra->SetLineWidth(2);
	f3_trace_ra->SetLineWidth(2);
	f4_trace_ra->SetLineWidth(2);
	f1_trace_ra->SetLineColor(1);
	f2_trace_ra->SetLineColor(2);
	f3_trace_ra->SetLineColor(3);
	f4_trace_ra->SetLineColor(4);

	new TCanvas();
	f1_trace_ra->Draw("HIST SAME");
	f2_trace_ra->Draw("HIST SAME");
	f3_trace_ra->Draw("HIST SAME");
	f4_trace_ra->Draw("HIST SAME");




}
#include "../analysis_tools.C"

void plot_filtered_magnet(){

	
	//January comparison
	TFile* f1 = TFile::Open("analysis/analysis_EC_jan21_B3043_H24.root");
	TFile* f2 = TFile::Open("analysis/analysis_EC_jan20_B3619_H16.root");
	TFile* f3 = TFile::Open("analysis/analysis_EC_jan19_B4353_H23.root");
	TFile* f4 = TFile::Open("analysis/analysis_EC_jan19_B5175_H15.root");
	vector<double> strengths = {3043, 3619, 4353, 5175};

	TProfile* f1_trace_calibrated = (TProfile*)f1->Get("trace_calibrated");
	TProfile* f2_trace_calibrated = (TProfile*)f2->Get("trace_calibrated");
	TProfile* f3_trace_calibrated = (TProfile*)f3->Get("trace_calibrated");
	TProfile* f4_trace_calibrated = (TProfile*)f4->Get("trace_calibrated");


	TH1D* f1_kick1_calibrated = ((TProfile*)f1->Get("trace_kick1_calibrated"))->ProjectionX("f1_kick1_calibrated");
	TH1D* f2_kick1_calibrated = ((TProfile*)f2->Get("trace_kick1_calibrated"))->ProjectionX("f2_kick1_calibrated");
	TH1D* f3_kick1_calibrated = ((TProfile*)f3->Get("trace_kick1_calibrated"))->ProjectionX("f3_kick1_calibrated");
	TH1D* f4_kick1_calibrated = ((TProfile*)f4->Get("trace_kick1_calibrated"))->ProjectionX("f4_kick1_calibrated");
	f1_kick1_calibrated->SetLineColor(1);
	f2_kick1_calibrated->SetLineColor(2);
	f3_kick1_calibrated->SetLineColor(3);
	f4_kick1_calibrated->SetLineColor(4);
	f1_kick1_calibrated->SetLineWidth(2);
	f2_kick1_calibrated->SetLineWidth(2);
	f3_kick1_calibrated->SetLineWidth(2);
	f4_kick1_calibrated->SetLineWidth(2);

	double fft_start = 0.03;
	double fft_end = 2;

	TH1D* f1_kick1_calibrated_lowpass = lowpass(f1_kick1_calibrated,15,fft_start,fft_end);
	TH1D* f2_kick1_calibrated_lowpass = lowpass(f2_kick1_calibrated,15,fft_start,fft_end);
	TH1D* f3_kick1_calibrated_lowpass = lowpass(f3_kick1_calibrated,15,fft_start,fft_end);
	TH1D* f4_kick1_calibrated_lowpass = lowpass(f4_kick1_calibrated,15,fft_start,fft_end);
	f1_kick1_calibrated_lowpass->SetLineColor(1);
	f2_kick1_calibrated_lowpass->SetLineColor(2);
	f3_kick1_calibrated_lowpass->SetLineColor(3);
	f4_kick1_calibrated_lowpass->SetLineColor(4);
	f1_kick1_calibrated_lowpass->SetLineWidth(2);
	f2_kick1_calibrated_lowpass->SetLineWidth(2);
	f3_kick1_calibrated_lowpass->SetLineWidth(2);
	f4_kick1_calibrated_lowpass->SetLineWidth(2);


	TH1D* f1_kick1_calibrated_FFT = doFFT(f1_kick1_calibrated,fft_start,fft_end);
	TH1D* f2_kick1_calibrated_FFT = doFFT(f2_kick1_calibrated,fft_start,fft_end);
	TH1D* f3_kick1_calibrated_FFT = doFFT(f3_kick1_calibrated,fft_start,fft_end);
	TH1D* f4_kick1_calibrated_FFT = doFFT(f4_kick1_calibrated,fft_start,fft_end);
	f1_kick1_calibrated_FFT->SetLineColor(1);
	f2_kick1_calibrated_FFT->SetLineColor(2);
	f3_kick1_calibrated_FFT->SetLineColor(3);
	f4_kick1_calibrated_FFT->SetLineColor(4);
	f1_kick1_calibrated_FFT->SetLineWidth(2);
	f2_kick1_calibrated_FFT->SetLineWidth(2);
	f3_kick1_calibrated_FFT->SetLineWidth(2);
	f4_kick1_calibrated_FFT->SetLineWidth(2);

	TH1D* f1_kick1_calibrated_lowpass_FFT = doFFT(f1_kick1_calibrated_lowpass,fft_start,fft_end);
	TH1D* f2_kick1_calibrated_lowpass_FFT = doFFT(f2_kick1_calibrated_lowpass,fft_start,fft_end);
	TH1D* f3_kick1_calibrated_lowpass_FFT = doFFT(f3_kick1_calibrated_lowpass,fft_start,fft_end);
	TH1D* f4_kick1_calibrated_lowpass_FFT = doFFT(f4_kick1_calibrated_lowpass,fft_start,fft_end);
	f1_kick1_calibrated_lowpass_FFT->SetLineColor(1);
	f2_kick1_calibrated_lowpass_FFT->SetLineColor(2);
	f3_kick1_calibrated_lowpass_FFT->SetLineColor(3);
	f4_kick1_calibrated_lowpass_FFT->SetLineColor(4);
	f1_kick1_calibrated_lowpass_FFT->SetLineWidth(2);
	f2_kick1_calibrated_lowpass_FFT->SetLineWidth(2);
	f3_kick1_calibrated_lowpass_FFT->SetLineWidth(2);
	f4_kick1_calibrated_lowpass_FFT->SetLineWidth(2);

	gStyle->SetOptStat(0);


	new TCanvas();
	f1_kick1_calibrated->GetXaxis()->SetRangeUser(-1,2);
	f2_kick1_calibrated->GetXaxis()->SetRangeUser(-1,2);
	f3_kick1_calibrated->GetXaxis()->SetRangeUser(-1,2);
	f4_kick1_calibrated->GetXaxis()->SetRangeUser(-1,2);
	f1_kick1_calibrated->GetYaxis()->SetRangeUser(-150,150);
	f2_kick1_calibrated->GetYaxis()->SetRangeUser(-150,150);
	f3_kick1_calibrated->GetYaxis()->SetRangeUser(-150,150);
	f4_kick1_calibrated->GetYaxis()->SetRangeUser(-150,150);
	f1_kick1_calibrated->Draw("HIST SAME");
	f2_kick1_calibrated->Draw("HIST SAME");
	f3_kick1_calibrated->Draw("HIST SAME");
	f4_kick1_calibrated->Draw("HIST SAME");

	new TCanvas();
	f1_kick1_calibrated_lowpass->GetXaxis()->SetRangeUser(-1,2);
	f2_kick1_calibrated_lowpass->GetXaxis()->SetRangeUser(-1,2);
	f3_kick1_calibrated_lowpass->GetXaxis()->SetRangeUser(-1,2);
	f4_kick1_calibrated_lowpass->GetXaxis()->SetRangeUser(-1,2);
	f1_kick1_calibrated_lowpass->GetYaxis()->SetRangeUser(-150,150);
	f2_kick1_calibrated_lowpass->GetYaxis()->SetRangeUser(-150,150);
	f3_kick1_calibrated_lowpass->GetYaxis()->SetRangeUser(-150,150);
	f4_kick1_calibrated_lowpass->GetYaxis()->SetRangeUser(-150,150);
	f1_kick1_calibrated_lowpass->Draw("HIST SAME");
	f2_kick1_calibrated_lowpass->Draw("HIST SAME");
	f3_kick1_calibrated_lowpass->Draw("HIST SAME");
	f4_kick1_calibrated_lowpass->Draw("HIST SAME");


	new TCanvas();
	f1_kick1_calibrated_FFT->Draw("HIST SAME");
	f2_kick1_calibrated_FFT->Draw("HIST SAME");
	f3_kick1_calibrated_FFT->Draw("HIST SAME");
	f4_kick1_calibrated_FFT->Draw("HIST SAME");

	new TCanvas();
	f1_kick1_calibrated_lowpass_FFT->Draw("HIST SAME");
	f2_kick1_calibrated_lowpass_FFT->Draw("HIST SAME");
	f3_kick1_calibrated_lowpass_FFT->Draw("HIST SAME");
	f4_kick1_calibrated_lowpass_FFT->Draw("HIST SAME");

}
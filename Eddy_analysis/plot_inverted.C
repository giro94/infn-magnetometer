#include "../analysis_tools.C"

void plot_inverted(){




	TFile* f1 = TFile::Open("kick_analysis_output_SD_R0_eddy_dec15_H25_B100_k777.root");
	TFile* f2 = TFile::Open("kick_analysis_output_SD_R0_eddy_dec16_H25_B100_k777_inverted.root");

	TH1D* h1 = ((TProfile*)f1->Get("trace_SNR2_timealigned"))->ProjectionX("h1");
	TH1D* h2 = ((TProfile*)f2->Get("trace_SNR2_timealigned"))->ProjectionX("h2");

	//Normalize blumlein to 50 mV
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");

	f_blumlein->SetParameters(60.0,-0.3,-1000.0);
	TFitResultPtr fit_1 = h1->Fit("f_blumlein","NQS","",-0.4,-0.2);
	if (fit_1>=0){
		double peak = fit_1->Parameter(0);
		h1->Scale(50./peak);
	}
	f_blumlein->SetParameters(-60.0,-0.3,-1000.0);
	TFitResultPtr fit_2 = h2->Fit("f_blumlein","NQS","",-0.4,-0.2);
	if (fit_2>=0){
		double peak = fit_2->Parameter(0);
		h2->Scale(-50./peak);
	}

	TH1D* hsum = (TH1D*)h1->Clone("hsum");
	hsum->Add(h2);
	//hsum->Scale(0.5);

	TH1D* hdiff = (TH1D*)h1->Clone("hdiff");
	hdiff->Add(h2,-1);
	hdiff->Scale(0.5);


	new TCanvas();
	h1->SetLineColor(kBlue);
	h2->SetLineColor(kRed);
	//h1->GetXaxis()->SetRangeUser(-0.5,2);
	h1->GetYaxis()->SetRangeUser(-60,60);
	//h2->GetXaxis()->SetRangeUser(-0.5,2);
	h2->GetYaxis()->SetRangeUser(-60,60);
	h1->Draw("HIST");
	h2->Draw("HIST SAME");
	gPad->SetGridy();


	new TCanvas();
	hsum->SetLineColor(kViolet);
	//hsum->GetXaxis()->SetRangeUser(-0.5,2);
	hsum->GetYaxis()->SetRangeUser(-60,60);
	hsum->Draw("HIST");

	new TCanvas();
	hdiff->SetLineColor(kOrange);
	//hdiff->GetXaxis()->SetRangeUser(-0.5,2);
	hdiff->GetYaxis()->SetRangeUser(-60,60);
	hdiff->Draw("HIST");




	new TCanvas();
	h1->Draw("HIST");
	hdiff->Draw("HIST SAME");	





}
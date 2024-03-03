#include "../analysis_tools.C"

void plot_BField_ONOFF(){

	TFile* f0_1 = TFile::Open("analysis/analysis_SD_R0_eddy_oct14_H25_nofilter_B0.root");
	TFile* f0_2 = TFile::Open("analysis/analysis_SD_R0_eddy_oct23_H25_B0_k777.root");

	TFile* f1_1 = TFile::Open("analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f1_2 = TFile::Open("analysis/analysis_EC_jan28_B5173_H25_Q00.root");

	TFile* f2 = TFile::Open("analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");


	TH1D* h1_kick1_B0_1 = ((TProfile*)f0_1->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B0_2 = ((TProfile*)f0_2->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B1_1 = ((TProfile*)f1_1->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B1_2 = ((TProfile*)f1_2->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_k0 = ((TProfile*)f2->Get("trace_kick1"))->ProjectionX();
	h1_kick1_k0->Scale(-1);
	//h1_kick1_B1_2->Scale(-1);

	TH1D* h1_trace_k0 = ((TProfile*)f2->Get("trace"))->ProjectionX();


	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1_kick1_B0_1_ra = smoothing(h1_kick1_B0_1,"");
	TH1D* h1_kick1_B0_2_ra = smoothing(h1_kick1_B0_2,"");
	TH1D* h1_kick1_B1_1_ra = smoothing(h1_kick1_B1_1,"");
	TH1D* h1_kick1_B1_2_ra = smoothing(h1_kick1_B1_2,"");
	TH1D* h1_kick1_k0_ra = smoothing(h1_kick1_k0,"");




	TH1D* h1_kick1_B0 = (TH1D*)h1_kick1_B0_1->Clone("h1_kick1_B0");
	h1_kick1_B0->Add(h1_kick1_B0_2);
	h1_kick1_B0->Scale(0.5);


	TH1D* h1_kick1_B1 = (TH1D*)h1_kick1_B1_1->Clone("h1_kick1_B1");
	h1_kick1_B1->Add(h1_kick1_B1_2);
	h1_kick1_B1->Scale(0.5);


	TH1D* h1_kick1_B0_ra = smoothing(h1_kick1_B0,"");
	TH1D* h1_kick1_B1_ra = smoothing(h1_kick1_B1,"");


	TH1D* h1_kick1_B1_1_k0 = (TH1D*)h1_kick1_B1_1->Clone("h1_kick1_B1_1_k0");
	h1_kick1_B1_1_k0->Add(h1_kick1_k0,-1);
	TH1D* h1_kick1_B1_1_k0_ra = smoothing(h1_kick1_B1_1_k0,"");

	TH1D* h1_kick1_B1_2_k0 = (TH1D*)h1_kick1_B1_2->Clone("h1_kick1_B1_2_k0");
	h1_kick1_B1_2_k0->Add(h1_kick1_k0,-1);
	TH1D* h1_kick1_B1_2_k0_ra = smoothing(h1_kick1_B1_2_k0,"");

	TH1D* h1_kick1_B1_k0 = (TH1D*)h1_kick1_B1->Clone("h1_kick1_B1_k0");
	h1_kick1_B1_k0->Add(h1_kick1_k0,-1);
	TH1D* h1_kick1_B1_k0_ra = smoothing(h1_kick1_B1_k0,"");



	double xmin = -1;
	double xmax = 8;
	int nbins = 8928;

	TH1D* h1_kick1_B0_resampled = new TH1D("h1_kick1_B0_resampled","",nbins,xmin,xmax);
	TH1D* h1_kick1_B0_ra_resampled = new TH1D("h1_kick1_B0_ra_resampled","",nbins,xmin,xmax);
	for (int bx=1; bx<=h1_kick1_B0_resampled->GetNbinsX(); bx++){
		double bin_center = h1_kick1_B0_resampled->GetBinCenter(bx);
		double yval = h1_kick1_B0->Interpolate(bin_center);
		h1_kick1_B0_resampled->SetBinContent(bx,yval);
		double yval2 = h1_kick1_B0_ra->Interpolate(bin_center);
		h1_kick1_B0_ra_resampled->SetBinContent(bx,yval2);
	}

	TH1D* h1_kick1_B1_B0 = (TH1D*)h1_kick1_B1->Clone("h1_kick1_B1_B0");
	h1_kick1_B1_B0->Add(h1_kick1_B0_resampled,-1);

	TH1D* h1_kick1_B1_B0_ra = smoothing(h1_kick1_B1_B0,"");

	TH1D* h1_kick1_B1_B0_ra0 = (TH1D*)h1_kick1_B1_ra->Clone("h1_kick1_B1_B0_ra0");
	h1_kick1_B1_B0_ra0->Add(h1_kick1_B0_ra_resampled,-1);


	TH2D* h1_kick1_B1_B0_ra_spectrograph = getSpectrograph(h1_kick1_B1_B0_ra, 0, 8, 0.5);

	TH2D* h1_trace_k0_spectrograph = getSpectrograph(h1_trace_k0, 0, 500, 10);

//////////////////////// Style

	h1_kick1_B0_1->SetTitle("Kick 1 (0 A) (Oct 14)");
	h1_kick1_B0_2->SetTitle("Kick 1 (0 A) (Oct 23)");
	h1_kick1_B1_1->SetTitle("Kick 1 (5173 A) (Jan 26)");
	h1_kick1_B1_2->SetTitle("Kick 1 (5173 A) (Jan 28)");
	h1_kick1_B0_1->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_2->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_1->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_2->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_1->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_2->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_1->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_2->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_1->SetLineWidth(2);
	h1_kick1_B0_2->SetLineWidth(2);
	h1_kick1_B1_1->SetLineWidth(2);
	h1_kick1_B1_2->SetLineWidth(2);
	h1_kick1_B0_1->SetLineColor(kBlue);
	h1_kick1_B0_2->SetLineColor(kRed);
	h1_kick1_B1_1->SetLineColor(kBlue);
	h1_kick1_B1_2->SetLineColor(kRed);

	h1_kick1_B0_1_ra->SetTitle("Kick 1 (0 A) (Oct 14)");
	h1_kick1_B0_2_ra->SetTitle("Kick 1 (0 A) (Oct 23)");
	h1_kick1_B1_1_ra->SetTitle("Kick 1 (5173 A) (Jan 26)");
	h1_kick1_B1_2_ra->SetTitle("Kick 1 (5173 A) (Jan 28)");
	h1_kick1_B0_1_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_2_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_1_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_2_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_1_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_2_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_1_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_2_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_1_ra->SetLineWidth(2);
	h1_kick1_B0_2_ra->SetLineWidth(2);
	h1_kick1_B1_1_ra->SetLineWidth(2);
	h1_kick1_B1_2_ra->SetLineWidth(2);
	h1_kick1_B0_1_ra->SetLineColor(kBlue);
	h1_kick1_B0_2_ra->SetLineColor(kRed);
	h1_kick1_B1_1_ra->SetLineColor(kBlue);
	h1_kick1_B1_2_ra->SetLineColor(kRed);



	h1_kick1_k0_ra->SetTitle("Kick 1 (5173 A, QWP) (Jan 25)");
	h1_kick1_k0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_k0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_k0_ra->SetLineWidth(2);
	h1_kick1_k0_ra->SetLineColor(kBlack);

	h1_kick1_B0->SetTitle("Kick 1 (0 A) (Oct 14 + 23)");
	h1_kick1_B1->SetTitle("Kick 1 (5173 A) (Jan 26 + 28)");
	h1_kick1_B0->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0->SetLineWidth(2);
	h1_kick1_B1->SetLineWidth(2);
	h1_kick1_B0->SetLineColor(kBlue);
	h1_kick1_B1->SetLineColor(kRed);



	h1_kick1_B0_ra->SetTitle("Kick 1 (0 A) (Oct 14 + 23)");
	h1_kick1_B1_ra->SetTitle("Kick 1 (5173 A) (Jan 26 + 28)");
	h1_kick1_B0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_ra->SetLineWidth(2);
	h1_kick1_B1_ra->SetLineWidth(2);
	h1_kick1_B0_ra->SetLineColor(kBlue);
	h1_kick1_B1_ra->SetLineColor(kRed);


	h1_kick1_B1_B0->SetTitle("Kick 1 (B - noB)");
	h1_kick1_B1_B0_ra->SetTitle("Kick 1 (B - noB)");
	h1_kick1_B1_B0_ra0->SetTitle("Kick 1 (B - noB)");
	h1_kick1_B1_B0->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_B0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_B0_ra0->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_B0->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_B0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_B0_ra0->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_B0->SetLineWidth(2);
	h1_kick1_B1_B0_ra->SetLineWidth(2);
	h1_kick1_B1_B0_ra0->SetLineWidth(2);
	h1_kick1_B1_B0->SetLineColor(kBlue);
	h1_kick1_B1_B0_ra->SetLineColor(kRed);
	h1_kick1_B1_B0_ra0->SetLineColor(kViolet);




//////////////////////// Drawing


	gStyle->SetOptStat(0);

	new TCanvas();
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B0_2_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	h1_kick1_B1_2_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_1_ra->SetLineColor(1);
	h1_kick1_B0_2_ra->SetLineColor(2);
	h1_kick1_B1_1_ra->SetLineColor(3);
	h1_kick1_B1_2_ra->SetLineColor(4);
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B0_2_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_1_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_2_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas();
	h1_kick1_B0_1->SetLineColor(1);
	h1_kick1_B0_2->SetLineColor(2);
	h1_kick1_B1_1->SetLineColor(3);
	h1_kick1_B1_2->SetLineColor(4);
	h1_kick1_B0_1->DrawCopy("HIST");
	h1_kick1_B0_2->DrawCopy("HIST SAME");
	h1_kick1_B1_1->DrawCopy("HIST SAME");
	h1_kick1_B1_2->DrawCopy("HIST SAME");
	gPad->BuildLegend();




	new TCanvas();
	h1_kick1_B0->DrawCopy("HIST");
	h1_kick1_B1->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas();
	h1_kick1_B0_1->SetLineColor(1);
	h1_kick1_B0_2->SetLineColor(2);
	h1_kick1_B0->SetLineColor(3);
	h1_kick1_B0_1->DrawCopy("HIST");
	h1_kick1_B0_2->DrawCopy("HIST SAME");
	h1_kick1_B0->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_1->SetLineColor(1);
	h1_kick1_B1_2->SetLineColor(2);
	h1_kick1_B1->SetLineColor(3);
	h1_kick1_B1_1->DrawCopy("HIST");
	h1_kick1_B1_2->DrawCopy("HIST SAME");
	h1_kick1_B1->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_1_ra->SetLineColor(1);
	h1_kick1_B0_2_ra->SetLineColor(2);
	h1_kick1_B0_ra->SetLineColor(3);
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B0_2_ra->DrawCopy("HIST SAME");
	h1_kick1_B0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_1_ra->SetLineColor(1);
	h1_kick1_B1_2_ra->SetLineColor(2);
	h1_kick1_B1_ra->SetLineColor(3);
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	h1_kick1_B1_2_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_B0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_B0_ra0->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B0_ra_resampled->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B1_B0_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B1_1_ra->Scale(-1);
	h1_kick1_B1_1_ra->SetLineColor(kRed);
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B1_2_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	h1_kick1_B1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_2_ra->DrawCopy("HIST");
	h1_kick1_B1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B1_ra->DrawCopy("HIST");
	h1_kick1_B1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();
	
	new TCanvas();
	h1_kick1_B0_ra_resampled->DrawCopy("HIST");
	h1_kick1_B1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas("","",1200,800);
	h1_kick1_B0_1_ra->SetLineWidth(2);
	h1_kick1_B0_1_ra->SetLineColor(kBlack);
	h1_kick1_B1_k0_ra->SetLineWidth(2);
	h1_kick1_B1_k0_ra->SetLineColor(kRed);
	h1_kick1_B0_1_ra->GetXaxis()->SetRangeUser(0,1);
	h1_kick1_B0_1_ra->GetYaxis()->SetRangeUser(-30,20);
	h1_kick1_B0_1_ra->SetTitle("Bfield OFF (Oct)");
	h1_kick1_B1_k0_ra->SetTitle("Bfield ON (Jan)");
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B1_k0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend(0.25,0.15,0.55,0.3);
	gPad->SetGridy();
	TLine* l_30 = new TLine(0.03,-30,0.03,20);
	TLine* l_700 = new TLine(0.7,-30,0.7,20);
	l_30->SetLineStyle(kDashed);
	l_700->SetLineStyle(kDashed);
	l_30->Draw();
	l_700->Draw();


	gStyle->SetPalette(kRainBow);
	new TCanvas();
	h1_kick1_B1_B0_ra_spectrograph->SetContour(100);
	h1_kick1_B1_B0_ra_spectrograph->DrawCopy("colz");


	new TCanvas();
	h1_trace_k0_spectrograph->SetContour(100);
	h1_trace_k0_spectrograph->DrawCopy("colz");




	TCanvas* can_qwp = new TCanvas("","",1400,1200);
	can_qwp->Divide(3,3);
	h1_kick1_B1_1_ra->SetLineColor(kBlue);
	h1_kick1_B1_2_ra->SetLineColor(kBlue);
	h1_kick1_B1_ra->SetLineColor(kBlue);
	h1_kick1_k0_ra->SetLineColor(kBlack);
	h1_kick1_B1_1_k0_ra->SetLineColor(kRed);
	h1_kick1_B1_2_k0_ra->SetLineColor(kRed);
	h1_kick1_B1_k0_ra->SetLineColor(kRed);
	h1_kick1_B1_1_ra->SetLineWidth(2);
	h1_kick1_B1_2_ra->SetLineWidth(2);
	h1_kick1_B1_ra->SetLineWidth(2);
	h1_kick1_k0_ra->SetLineWidth(2);
	h1_kick1_B1_1_k0_ra->SetLineWidth(2);
	h1_kick1_B1_2_k0_ra->SetLineWidth(2);
	h1_kick1_B1_k0_ra->SetLineWidth(2);
	h1_kick1_B1_1_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_2_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_k0_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_1_k0_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_2_k0_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_k0_ra->GetXaxis()->SetRangeUser(1,1.5);
	h1_kick1_B1_1_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_B1_2_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_B1_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_k0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_B1_1_k0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_B1_2_k0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_B1_k0_ra->GetYaxis()->SetRangeUser(-50,50);
	can_qwp->cd(1);
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(2);
	h1_kick1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(3);
	h1_kick1_B1_1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(4);
	h1_kick1_B1_2_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(5);
	h1_kick1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(6);
	h1_kick1_B1_2_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(7);
	h1_kick1_B1_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(8);
	h1_kick1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp->cd(9);
	h1_kick1_B1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();


	TCanvas* can_qwp2 = new TCanvas("","",1400,1200);
	can_qwp2->Divide(2,3);
	can_qwp2->cd(1);
	h1_kick1_B1_1_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->SetGridy();
	can_qwp2->cd(2);
	h1_kick1_B1_1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp2->cd(3);
	h1_kick1_B1_2_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->SetGridy();
	can_qwp2->cd(4);
	h1_kick1_B1_2_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();
	can_qwp2->cd(5);
	h1_kick1_B1_ra->DrawCopy("HIST");
	h1_kick1_k0_ra->DrawCopy("HIST SAME");
	gPad->SetGridy();
	can_qwp2->cd(6);
	h1_kick1_B1_k0_ra->DrawCopy("HIST");
	gPad->SetGridy();

}
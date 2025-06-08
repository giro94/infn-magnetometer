#include "../analysis_tools.C"

void plot_kicks_paper(){

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10;

	TFile* fin = TFile::Open("INFN_UMass.root");

	TH1D* h1_kick1_R0 = (TH1D*)fin->Get("h1_kick1_R0");
	TH1D* h1_kick1_R0_ra = (TH1D*)fin->Get("h1_kick1_R0_ra");
	TH1D* h1_kick1_R1 = (TH1D*)fin->Get("h1_kick1_R1");
	TH1D* h1_kick1_R1_ra = (TH1D*)fin->Get("h1_kick1_R1_ra");


	TH1D* h1_K1_R0 = (TH1D*)fin->Get("h1_K1_R0");
	TH1D* h1_K3_R0 = (TH1D*)fin->Get("h1_K3_R0");
	TH1D* h1_K1_R0_2022 = (TH1D*)fin->Get("h1_K1_R0_2022");
	TH1D* h1_K1_R0_ra = smoothing(h1_K1_R0,"");
	TH1D* h1_K3_R0_ra = smoothing(h1_K3_R0,"");
	TH1D* h1_K1_R0_2022_ra = smoothing(h1_K1_R0_2022,"");

	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");
	f_blumlein->SetParameters(125.0,-0.3,-1000.0);
	double xmin_blum = -0.42;
	double xmax_blum = -0.22;

	/////////Normalize INFN according to final numbers
	f_blumlein->SetParameters(125.0,-0.32,-3000.0);
	double blum_R0 = 126.397;
	double blum_R1 = 158.218;
	TFitResultPtr fit_blum_R0 = h1_kick1_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_R1 = h1_kick1_R1_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	double blum_fit_R0 = fit_blum_R0->Parameter(0);
	double blum_fit_R1 = fit_blum_R1->Parameter(0);
	h1_kick1_R0_ra->Scale(blum_R0/blum_fit_R0);
	h1_kick1_R1_ra->Scale(blum_R1/blum_fit_R1);
	/////////

	/////////Normalize UMass according to 30/01/2025 email from David
	f_blumlein->SetParameters(125.0,-0.32,-3000.0);
	double blum_K1_R0 = 123.0;
	double blum_K3_R0 = 126.3;
	double blum_K1_R0_2022 = 124.8;
	TFitResultPtr fit_blum_K1 = h1_K1_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K3 = h1_K3_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K1_2022 = h1_K1_R0_2022_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	double blum_fit_K1 = fit_blum_K1->Parameter(0);
	double blum_fit_K3 = fit_blum_K3->Parameter(0);
	double blum_fit_K1_2022 = fit_blum_K1_2022->Parameter(0);
	h1_K1_R0_ra->Scale(blum_K1_R0/blum_fit_K1);
	h1_K3_R0_ra->Scale(blum_K3_R0/blum_fit_K3);
	h1_K1_R0_2022_ra->Scale(blum_K1_R0_2022/blum_fit_K1_2022);
	/////////

	h1_K1_R0_ra->SetTitle(Form("%s (smoothed)", h1_K1_R0->GetTitle()));
	h1_K3_R0_ra->SetTitle(Form("%s (smoothed)", h1_K3_R0->GetTitle()));
	h1_K1_R0_2022_ra->SetTitle(Form("%s (smoothed)", h1_K1_R0_2022->GetTitle()));

	double f_kickers = (53.1+53.0+55.0)/(3*55.0);

	double k1_to_k3 = 55.0/53.1;

	TH1D* h1_K1_R0_norm = (TH1D*)h1_K1_R0->Clone("h1_K1_R0_norm");
	h1_K1_R0_norm->Scale(k1_to_k3);
	h1_K1_R0_norm->SetTitle(Form("UMass %s", h1_K1_R0->GetTitle()));

	TH1D* h1_K1_R0_2022_norm = (TH1D*)h1_K1_R0_2022->Clone("h1_K1_R0_2022_norm");
	h1_K1_R0_2022_norm->Scale(k1_to_k3);
	h1_K1_R0_2022_norm->SetTitle(Form("UMass %s", h1_K1_R0_2022->GetTitle()));

	TH1D* h1_K1_R0_ra_norm = (TH1D*)h1_K1_R0_ra->Clone("h1_K1_R0_ra_norm");
	h1_K1_R0_ra_norm->Scale(k1_to_k3);
	h1_K1_R0_ra_norm->SetTitle(Form("UMass %s", h1_K1_R0->GetTitle()));

	TH1D* h1_K1_R0_2022_ra_norm = (TH1D*)h1_K1_R0_2022_ra->Clone("h1_K1_R0_2022_ra_norm");
	h1_K1_R0_2022_ra_norm->Scale(k1_to_k3);
	h1_K1_R0_2022_ra_norm->SetTitle(Form("UMass %s", h1_K1_R0_2022->GetTitle()));


	TH1D* h1_kick1_R1_ra_norm = (TH1D*)h1_kick1_R1_ra->Clone("h1_kick1_R1_ra_norm");
	h1_kick1_R1_ra_norm->Scale(h1_kick1_R0_ra->Interpolate(0.03)/h1_kick1_R1_ra_norm->Interpolate(0.03));
	h1_kick1_R1_ra_norm->SetTitle("INFN R1 normalized to R0");


	//Resample UMass K1 & K3 to match INFN binning (slightly larger)
	TH1D* h1_K1_R0_resampled = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R0_resampled");
	TH1D* h1_K3_R0_resampled = (TH1D*)h1_kick1_R0_ra->Clone("h1_K3_R0_resampled");
	TH1D* h1_K1_R0_2022_resampled = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R0_2022_resampled");
	h1_K1_R0_resampled->Reset();
	h1_K3_R0_resampled->Reset();
	h1_K1_R0_2022_resampled->Reset();
	for (int bn=1; bn<=h1_kick1_R0_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		double y1 = h1_K1_R0_ra_norm->Interpolate(x);
		h1_K1_R0_resampled->SetBinContent(bn,y1);
		double y3 = h1_K3_R0_ra->Interpolate(x);
		h1_K3_R0_resampled->SetBinContent(bn,y3);
		double y2 = h1_K1_R0_2022_ra_norm->Interpolate(x);
		h1_K1_R0_2022_resampled->SetBinContent(bn,y2);
	}

	gStyle->SetOptStat(0);

	new TCanvas("","smoothed, rescaled, normalized",1200,1000);
	h1_K1_R0_resampled->SetTitle("UMass K1 (rescaled for K3)");
	h1_K3_R0_resampled->SetTitle("UMass K3");
	h1_K1_R0_2022_resampled->SetTitle("UMass K1 2022 (rescaled for K3)");
	h1_kick1_R0_ra->SetTitle("INFN K3");
	h1_kick1_R1_ra_norm->SetTitle("INFN K3 R1 (rescaled for R0)");

	h1_K1_R0_resampled->SetLineWidth(2);
	h1_K3_R0_resampled->SetLineWidth(2);
	h1_K1_R0_2022_resampled->SetLineWidth(2);
	h1_kick1_R0_ra->SetLineWidth(2);
	h1_kick1_R1_ra_norm->SetLineWidth(2);

	h1_K1_R0_resampled->SetLineColor(kBlack);
	h1_K3_R0_resampled->SetLineColor(kBlue);
	h1_K1_R0_2022_resampled->SetLineColor(kGreen);
	h1_kick1_R0_ra->SetLineColor(kRed);
	h1_kick1_R1_ra_norm->SetLineColor(kViolet);

	double ylow = -25;
	double yhigh = 15;

	h1_K1_R0_resampled->GetXaxis()->SetRangeUser(0,0.8);
	h1_K1_R0_resampled->GetYaxis()->SetRangeUser(ylow,yhigh);

	h1_K1_R0_resampled->Draw("HIST");
	h1_K3_R0_resampled->Draw("HIST SAME");
	h1_K1_R0_2022_resampled->Draw("HIST SAME");
	h1_kick1_R0_ra->Draw("HIST SAME");
	//h1_kick1_R1_ra_norm->Draw("HIST SAME");

	TGraph* g_band = new TGraph();
	for (int bn=1; bn<=h1_kick1_R0_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		double y1 = h1_K1_R0_resampled->GetBinContent(bn);
		double y2 = h1_K3_R0_resampled->GetBinContent(bn);
		double y3 = h1_kick1_R0_ra->GetBinContent(bn);
		double y4 = h1_K1_R0_2022_resampled->GetBinContent(bn);

		double ymax = max({y1,y2,y3,y4});
		g_band->AddPoint(x,ymax);
	}
	//Return to 0
	for (int bn=h1_kick1_R0_ra->GetNbinsX(); bn>=1; bn--){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		double y1 = h1_K1_R0_resampled->GetBinContent(bn);
		double y2 = h1_K3_R0_resampled->GetBinContent(bn);
		double y3 = h1_kick1_R0_ra->GetBinContent(bn);
		double y4 = h1_K1_R0_2022_resampled->GetBinContent(bn);

		double ymin = min({y1,y2,y3,y4});
		g_band->AddPoint(x,ymin);
	}

	g_band->SetTitle("Band");
	g_band->SetLineWidth(1);
	g_band->SetLineColor(kRed);
	g_band->SetFillColor(kRed);
	g_band->SetFillStyle(3003);
	g_band->Draw("F");

	gPad->SetGridy();
	gPad->BuildLegend(0.35,0.25,0.65,0.45);

	TLine* l30 = new TLine(0.03,ylow,0.03,yhigh);
	TLine* l700 = new TLine(0.7,ylow,0.7,yhigh);
	l30->SetLineWidth(2);
	l700->SetLineWidth(2);
	l30->SetLineStyle(kDashed);
	l700->SetLineStyle(kDashed);
	l30->Draw("SAME");
	l700->Draw("SAME");


}
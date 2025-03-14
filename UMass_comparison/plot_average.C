#include "../analysis_tools.C"

void plot_average(){

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10;

	TFile* fin = TFile::Open("INFN_UMass.root");


	TH1D* h1_kick1_R0 = (TH1D*)fin->Get("h1_kick1_R0");
	TH1D* h1_kick1_R0_ra = (TH1D*)fin->Get("h1_kick1_R0_ra");
	TH1D* h1_kick1_R1 = (TH1D*)fin->Get("h1_kick1_R1");
	TH1D* h1_kick1_R1_ra = (TH1D*)fin->Get("h1_kick1_R1_ra");


	TH1D* h1_K1_R0 = (TH1D*)fin->Get("h1_K1_R0");
	TH1D* h1_K3_R0 = (TH1D*)fin->Get("h1_K3_R0");
	TH1D* h1_K1_R0_ra = smoothing(h1_K1_R0,"");
	TH1D* h1_K3_R0_ra = smoothing(h1_K3_R0,"");

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
	TFitResultPtr fit_blum_K1 = h1_K1_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K3 = h1_K3_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	double blum_fit_K1 = fit_blum_K1->Parameter(0);
	double blum_fit_K3 = fit_blum_K3->Parameter(0);
	h1_K1_R0_ra->Scale(blum_K1_R0/blum_fit_K1);
	h1_K3_R0_ra->Scale(blum_K3_R0/blum_fit_K3);
	/////////

	h1_K1_R0_ra->SetTitle(Form("%s (smoothed)", h1_K1_R0->GetTitle()));
	h1_K3_R0_ra->SetTitle(Form("%s (smoothed)", h1_K3_R0->GetTitle()));

	double f_kickers = (53.1+53.0+55.0)/(3*55.0);

	double k1_to_k3 = 55.0/53.1;

	TH1D* h1_K1_R0_norm = (TH1D*)h1_K1_R0->Clone("h1_K1_R0_norm");
	h1_K1_R0_norm->Scale(k1_to_k3);
	h1_K1_R0_norm->SetTitle(Form("UMass %s", h1_K1_R0->GetTitle()));

	TH1D* h1_K1_R0_ra_norm = (TH1D*)h1_K1_R0_ra->Clone("h1_K1_R0_ra_norm");
	h1_K1_R0_ra_norm->Scale(k1_to_k3);
	h1_K1_R0_ra_norm->SetTitle(Form("UMass %s", h1_K1_R0->GetTitle()));


	TH1D* h1_kick1_R1_ra_norm = (TH1D*)h1_kick1_R1_ra->Clone("h1_kick1_R1_ra_norm");
	h1_kick1_R1_ra_norm->Scale(h1_kick1_R0_ra->Interpolate(0.03)/h1_kick1_R1_ra_norm->Interpolate(0.03));
	h1_kick1_R1_ra_norm->SetTitle("INFN R1 normalized to R0");


	//Resample UMass K1 & K3 to match INFN binning (slightly larger)
	TH1D* h1_K1_R0_resampled = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R0_resampled");
	TH1D* h1_K3_R0_resampled = (TH1D*)h1_kick1_R0_ra->Clone("h1_K3_R0_resampled");
	h1_K1_R0_resampled->Reset();
	h1_K3_R0_resampled->Reset();
	for (int bn=1; bn<=h1_kick1_R0_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		double y1 = h1_K1_R0_ra_norm->Interpolate(x);
		h1_K1_R0_resampled->SetBinContent(bn,y1);
		double y3 = h1_K3_R0_ra->Interpolate(x);
		h1_K3_R0_resampled->SetBinContent(bn,y3);
	}



	new TCanvas("","Original");
	h1_K1_R0->SetTitle("UMass K1");
	h1_K3_R0->SetTitle("UMass K3");
	h1_kick1_R0_ra->SetTitle("INFN K3");

	h1_K1_R0->SetLineWidth(2);
	h1_K3_R0->SetLineWidth(2);
	h1_kick1_R0_ra->SetLineWidth(2);

	h1_K1_R0->SetLineColor(kBlack);
	h1_K3_R0->SetLineColor(kBlue);
	h1_kick1_R0_ra->SetLineColor(kRed);

	h1_K1_R0->Draw("HIST");
	h1_K3_R0->Draw("HIST SAME");
	h1_kick1_R0_ra->Draw("HIST SAME");


	gStyle->SetOptStat(0);

	new TCanvas("","smoothed, rescaled, normalized",1200,1000);
	h1_K1_R0_resampled->SetTitle("UMass K1 (rescaled for K3)");
	h1_K3_R0_resampled->SetTitle("UMass K3");
	h1_kick1_R0_ra->SetTitle("INFN K3");
	h1_kick1_R1_ra_norm->SetTitle("INFN K3 R1 (rescaled for R0)");

	h1_K1_R0_resampled->SetLineWidth(2);
	h1_K3_R0_resampled->SetLineWidth(2);
	h1_kick1_R0_ra->SetLineWidth(2);
	h1_kick1_R1_ra_norm->SetLineWidth(2);

	h1_K1_R0_resampled->SetLineColor(kBlack);
	h1_K3_R0_resampled->SetLineColor(kBlue);
	h1_kick1_R0_ra->SetLineColor(kRed);
	h1_kick1_R1_ra_norm->SetLineColor(kViolet);

	double ylow = -25;
	double yhigh = 15;

	h1_K1_R0_resampled->GetXaxis()->SetRangeUser(0,0.8);
	h1_K1_R0_resampled->GetYaxis()->SetRangeUser(ylow,yhigh);

	h1_K1_R0_resampled->Draw("HIST");
	h1_K3_R0_resampled->Draw("HIST SAME");
	h1_kick1_R0_ra->Draw("HIST SAME");
	h1_kick1_R1_ra_norm->Draw("HIST SAME");

	TH1D* h1_lower = (TH1D*)h1_kick1_R0_ra->Clone("h1_lower");
	TH1D* h1_middle = (TH1D*)h1_kick1_R0_ra->Clone("h1_middle");
	TH1D* h1_upper = (TH1D*)h1_kick1_R0_ra->Clone("h1_upper");
	h1_lower->Reset();
	h1_middle->Reset();
	h1_upper->Reset();
	TGraph* g_band = new TGraph();
	for (int bn=1; bn<=h1_kick1_R0_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		//if (x<0.03 || x>0.7) continue;

		double y1 = h1_K1_R0_resampled->GetBinContent(bn);
		double y2 = h1_K3_R0_resampled->GetBinContent(bn);
		double y3 = h1_kick1_R0_ra->GetBinContent(bn);

		double ymin = min({y1,y2,y3});
		double ymax = max({y1,y2,y3});
		double ymiddle = 0.5*(ymin+ymax);
		h1_lower->SetBinContent(bn,ymin);
		h1_middle->SetBinContent(bn,ymiddle);
		h1_upper->SetBinContent(bn,ymax);

		g_band->AddPoint(x,ymax);
	}
	//Return to 0 for the band
	for (int bn=h1_kick1_R0_ra->GetNbinsX(); bn>=1; bn--){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		//if (x<0.03 || x>0.7) continue;

		double y1 = h1_K1_R0_resampled->GetBinContent(bn);
		double y2 = h1_K3_R0_resampled->GetBinContent(bn);
		double y3 = h1_kick1_R0_ra->GetBinContent(bn);

		double ymin = min({y1,y2,y3});
		g_band->AddPoint(x,ymin);
	}

	g_band->SetTitle("Band");
	g_band->SetLineWidth(1);
	g_band->SetLineColor(kRed);
	g_band->SetFillColor(kRed);
	g_band->SetFillStyle(3003);
	g_band->Draw("F");

	//K1 center: Blumlein = 123.0 +/- 1.3 mG.
	//K3 center: Blumlein = 126.3 +/- 1.3 mG.
	//INFN R0: Blumlein = 126.4 +/- 2.9 mG.
	double er1 = (1.3/123.0);
	double er2 = (1.3/126.3);
	double er3 = (2.9/126.4);
	cout<<er1<<" "<<er2<<" "<<er3<<"\n";
	double w1 = 1./(er1*er1);
	double w2 = 1./(er2*er2);
	double w3 = 1./(er3*er3);
	double wtot = w1+w2+w3;
	w1 /= wtot;
	w2 /= wtot;
	w3 /= wtot;
	cout<<w1<<" "<<w2<<" "<<w3<<"\n";

	TH1D* h1_average = (TH1D*)h1_kick1_R0_ra->Clone("h1_average");
	TH1D* h1_average_w = (TH1D*)h1_kick1_R0_ra->Clone("h1_average_w");
	h1_average->Reset();
	h1_average_w->Reset();
	for (int bn=1; bn<=h1_average->GetNbinsX(); bn++){
		double x = h1_average->GetBinCenter(bn);

		double y1 = h1_K1_R0_resampled->GetBinContent(bn);
		double y2 = h1_K3_R0_resampled->GetBinContent(bn);
		double y3 = h1_kick1_R0_ra->GetBinContent(bn);

		double avg = (y1+y2+y3)/3.;
		double avg_w = ((y1*w1) + (y2*w2) + (y3*w3))/(w1+w2+w3);

		h1_average->SetBinContent(bn,avg);
		h1_average_w->SetBinContent(bn,avg_w);
	}

	h1_average->SetTitle("Arithmetic average");
	h1_average->SetLineWidth(2);
	h1_average->SetLineColor(kGreen);

	h1_average_w->SetTitle("Weighted average");
	h1_average_w->SetLineWidth(2);
	h1_average_w->SetLineColor(kGreen);

	//h1_average->Draw("HIST SAME");

	h1_lower->SetTitle("Lower transient");
	h1_lower->SetLineWidth(2);
	h1_lower->SetLineColor(kRed);

	h1_upper->SetTitle("Upper transient");
	h1_upper->SetLineWidth(2);
	h1_upper->SetLineColor(kRed);

	h1_middle->SetTitle("Middle transient");
	h1_middle->SetLineWidth(2);
	h1_middle->SetLineColor(kGreen);
	//h1_middle->Draw("HIST SAME");

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

	new TCanvas("","lower, middle, upper",1200,1000);
	h1_lower->GetXaxis()->SetRangeUser(0,0.8);
	h1_lower->GetYaxis()->SetRangeUser(ylow,yhigh);
	h1_lower->Draw("HIST");
	h1_middle->Draw("HIST SAME");
	h1_upper->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.35,0.25,0.65,0.45);
	l30->Draw("SAME");
	l700->Draw("SAME");


	new TCanvas("","lower, middle, upper",1200,1000);
	h1_average->GetXaxis()->SetRangeUser(0,0.8);
	h1_average->GetYaxis()->SetRangeUser(ylow,yhigh);
	h1_average->SetLineColor(kGreen);
	h1_average_w->SetLineColor(kViolet);
	h1_average->Draw("HIST");
	h1_average_w->Draw("HIST SAME");
	g_band->Draw("LF");
	gPad->SetGridy();
	gPad->BuildLegend(0.35,0.25,0.65,0.45);
	l30->Draw("SAME");
	l700->Draw("SAME");


	TFile* fout = new TFile("INFN_UMass_average.root","recreate");
	h1_K1_R0_resampled->Write("h1_UMass_K1_R0");
	h1_K3_R0_resampled->Write("h1_UMass_K3_R0");
	h1_kick1_R0_ra->Write("h1_INFN_R0");
	h1_kick1_R1_ra_norm->Write("h1_INFN_R1_R0norm");
	h1_lower->Write("h1_lower");
	h1_middle->Write("h1_middle");
	h1_upper->Write("h1_upper");
	h1_average->Write("h1_average");
	h1_average_w->Write("h1_average_w");
	fout->Write();
	fout->Close();

}
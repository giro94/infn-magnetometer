#include "../analysis_tools.C"

void plot_INFN_UMass(){

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10;

	TFile* fin = TFile::Open("INFN_UMass.root");


	TH1D* h1_kick1_R0_p = (TH1D*)fin->Get("h1_kick1_R0_p");
	TH1D* h1_kick1_R0_p_ra = (TH1D*)fin->Get("h1_kick1_R0_p_ra");
	TH1D* h1_kick1_R0 = (TH1D*)fin->Get("h1_kick1_R0");
	TH1D* h1_kick1_R0_ra = (TH1D*)fin->Get("h1_kick1_R0_ra");
	TH1D* h1_kick1_R1 = (TH1D*)fin->Get("h1_kick1_R1");
	TH1D* h1_kick1_R1_ra = (TH1D*)fin->Get("h1_kick1_R1_ra");


	TH1D* h1_K1_R0 = (TH1D*)fin->Get("h1_K1_R0");
	TH1D* h1_K1_R3p2 = (TH1D*)fin->Get("h1_K1_R3p2");
	TH1D* h1_K1_R6p6 = (TH1D*)fin->Get("h1_K1_R6p6");
	TH1D* h1_K1_Rm6p6 = (TH1D*)fin->Get("h1_K1_Rm6p6");
	TH1D* h1_K3_R0 = (TH1D*)fin->Get("h1_K3_R0");

	TGraphErrors* g_INFN_blumlein = new TGraphErrors();
	TGraphErrors* g_INFN_amp30 = new TGraphErrors();
	TGraphErrors* g_UMass_blumlein = new TGraphErrors();
	TGraphErrors* g_UMass_amp30 = new TGraphErrors();

	TH1D* h1_K1_R0_ra = smoothing(h1_K1_R0,"");
	TH1D* h1_K1_R3p2_ra = smoothing(h1_K1_R3p2,"");
	TH1D* h1_K1_R6p6_ra = smoothing(h1_K1_R6p6,"");
	TH1D* h1_K1_Rm6p6_ra = smoothing(h1_K1_Rm6p6,"");
	TH1D* h1_K3_R0_ra = smoothing(h1_K3_R0,"");

	h1_K1_R0_ra->SetTitle(Form("%s (smoothed)", h1_K1_R0->GetTitle()));
	h1_K1_R3p2_ra->SetTitle(Form("%s (smoothed)", h1_K1_R3p2->GetTitle()));
	h1_K1_R6p6_ra->SetTitle(Form("%s (smoothed)", h1_K1_R6p6->GetTitle()));
	h1_K1_Rm6p6_ra->SetTitle(Form("%s (smoothed)", h1_K1_Rm6p6->GetTitle()));
	h1_K3_R0_ra->SetTitle(Form("%s (smoothed)", h1_K3_R0->GetTitle()));

	TH1D* h1_K1_R0_ra_norm = (TH1D*)h1_K1_R0_ra->Clone("h1_K1_R0_ra_norm");
	TH1D* h1_K1_R3p2_ra_norm = (TH1D*)h1_K1_R3p2_ra->Clone("h1_K1_R3p2_ra_norm");
	TH1D* h1_K1_R6p6_ra_norm = (TH1D*)h1_K1_R6p6_ra->Clone("h1_K1_R6p6_ra_norm");
	TH1D* h1_K1_Rm6p6_ra_norm = (TH1D*)h1_K1_Rm6p6_ra->Clone("h1_K1_Rm6p6_ra_norm");
	TH1D* h1_K3_R0_ra_norm = (TH1D*)h1_K3_R0_ra->Clone("h1_K3_R0_ra_norm");

	TH1D* h1_kick1_R0_ra_norm = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_norm");
	TH1D* h1_kick1_R1_ra_norm = (TH1D*)h1_kick1_R1_ra->Clone("h1_kick1_R1_ra_norm");

	h1_K1_R0_ra_norm->Scale(1./h1_K1_R0_ra->Interpolate(-0.32));
	h1_K1_R3p2_ra_norm->Scale(1./h1_K1_R0_ra->Interpolate(-0.32));
	h1_K1_R6p6_ra_norm->Scale(1./h1_K1_R0_ra->Interpolate(-0.32));
	h1_K1_Rm6p6_ra_norm->Scale(1./h1_K1_R0_ra->Interpolate(-0.32));
	h1_K3_R0_ra_norm->Scale(1./h1_K3_R0_ra->Interpolate(-0.32));
	h1_kick1_R0_ra_norm->Scale(1./h1_kick1_R0_ra->Interpolate(-0.32));
	h1_kick1_R1_ra_norm->Scale(1./h1_kick1_R0_ra->Interpolate(-0.32));


	h1_K1_R0_ra_norm->SetTitle(Form("UMass %s", h1_K1_R0->GetTitle()));
	h1_K1_R3p2_ra_norm->SetTitle(Form("UMass %s", h1_K1_R3p2->GetTitle()));
	h1_K1_R6p6_ra_norm->SetTitle(Form("UMass %s", h1_K1_R6p6->GetTitle()));
	h1_K1_Rm6p6_ra_norm->SetTitle(Form("UMass %s", h1_K1_Rm6p6->GetTitle()));
	h1_K3_R0_ra_norm->SetTitle(Form("UMass %s", h1_K3_R0->GetTitle()));







	g_INFN_blumlein->AddPoint(0,h1_kick1_R0_ra->Interpolate(-0.32)/h1_kick1_R0_ra->Interpolate(-0.32));
	g_INFN_blumlein->AddPoint(17.5,h1_kick1_R1_ra->Interpolate(-0.32)/h1_kick1_R0_ra->Interpolate(-0.32));

	//Use averages!
	double blum_R0 = 126.397;
	double blum_R0err = 2.87448;
	double blum_R0err_relative = blum_R0err/blum_R0;
	double blum_R1 = 158.218;
	double blum_R1err = 3.59815;
	double blum_R1err_relative = blum_R1err/blum_R1;
	double blum_ratio = blum_R1/blum_R0;
	double blum_ratio_error = blum_ratio*sqrt(blum_R0err_relative*blum_R0err_relative + blum_R1err_relative*blum_R1err_relative);
	g_INFN_blumlein->SetPoint(0,0,1);
	g_INFN_blumlein->SetPointError(0,2,blum_R0err_relative);
	g_INFN_blumlein->SetPoint(1,17.5,blum_ratio);
	g_INFN_blumlein->SetPointError(1,2,blum_ratio_error);
	cout<<"blumlein ratio: "<<blum_ratio<<" +- "<<blum_ratio_error<<"\n";

	cout<<"INFN blum at K3, R0: "<<h1_kick1_R0_ra->Interpolate(-0.32)<<" mG\n";
	cout<<"UMass blum at K1, R0: "<<h1_K1_R0_ra->Interpolate(-0.32)<<" mG\n";
	g_UMass_blumlein->AddPoint(-6.6,h1_K1_Rm6p6_ra->Interpolate(-0.32)/h1_K1_R0_ra->Interpolate(-0.32));
	g_UMass_blumlein->AddPoint(0,h1_K1_R0_ra->Interpolate(-0.32)/h1_K1_R0_ra->Interpolate(-0.32));
	g_UMass_blumlein->AddPoint(3.2,h1_K1_R3p2_ra->Interpolate(-0.32)/h1_K1_R0_ra->Interpolate(-0.32));
	g_UMass_blumlein->AddPoint(6.6,h1_K1_R6p6_ra->Interpolate(-0.32)/h1_K1_R0_ra->Interpolate(-0.32));

	g_INFN_blumlein->SetMarkerStyle(20);
	g_UMass_blumlein->SetMarkerStyle(20);
	g_INFN_blumlein->SetMarkerColor(kBlue);
	g_UMass_blumlein->SetMarkerColor(kRed);
	g_INFN_blumlein->SetMarkerSize(1.5);
	g_UMass_blumlein->SetMarkerSize(1.5);


	g_INFN_amp30->AddPoint(0,h1_kick1_R0_ra->Interpolate(0.03)/h1_kick1_R0_ra->Interpolate(0.03));
	g_INFN_amp30->AddPoint(17.5,h1_kick1_R1_ra->Interpolate(0.03)/h1_kick1_R0_ra->Interpolate(0.03));

	//Use averages!
	double trans_R0 = 16.2271;
	double trans_R0err = 0.482309;
	double trans_R0err_relative = trans_R0err/trans_R0;
	double trans_R1 = 35.5134;
	double trans_R1err = 1.05555;
	double trans_R1err_relative = trans_R1err/trans_R1;
	double trans_ratio = trans_R1/trans_R0;
	double trans_ratio_error = trans_ratio*sqrt(trans_R0err_relative*trans_R0err_relative + blum_R1err_relative*blum_R1err_relative);
	g_INFN_amp30->SetPoint(0,0,1);
	g_INFN_amp30->SetPointError(0,2,trans_R0err_relative);
	g_INFN_amp30->SetPoint(1,17.5,trans_ratio);
	g_INFN_amp30->SetPointError(1,2,trans_ratio_error);
	cout<<"transient ratio: "<<trans_ratio<<" +- "<<trans_ratio_error<<"\n";

	g_UMass_amp30->AddPoint(-6.6,h1_K1_Rm6p6_ra->Interpolate(0.03)/h1_K1_R0_ra->Interpolate(0.03));
	g_UMass_amp30->AddPoint(0,h1_K1_R0_ra->Interpolate(0.03)/h1_K1_R0_ra->Interpolate(0.03));
	g_UMass_amp30->AddPoint(3.2,h1_K1_R3p2_ra->Interpolate(0.03)/h1_K1_R0_ra->Interpolate(0.03));
	g_UMass_amp30->AddPoint(6.6,h1_K1_R6p6_ra->Interpolate(0.03)/h1_K1_R0_ra->Interpolate(0.03));

	g_INFN_amp30->SetMarkerStyle(22);
	g_UMass_amp30->SetMarkerStyle(22);
	g_INFN_amp30->SetMarkerColor(kBlue);
	g_UMass_amp30->SetMarkerColor(kRed);
	g_INFN_amp30->SetMarkerSize(1.5);
	g_UMass_amp30->SetMarkerSize(1.5);


	g_INFN_blumlein->SetTitle("INFN (blumlein)");
	g_UMass_blumlein->SetTitle("UMass (blumlein)");
	g_INFN_amp30->SetTitle("INFN (transient at 30 #mus)");
	g_UMass_amp30->SetTitle("UMass (transient at 30 #mus)");

	g_INFN_blumlein->GetXaxis()->SetLimits(-25,25);
	g_UMass_blumlein->GetXaxis()->SetLimits(-25,25);
	g_INFN_amp30->GetXaxis()->SetLimits(-25,25);
	g_UMass_amp30->GetXaxis()->SetLimits(-25,25);
	g_INFN_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_UMass_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_INFN_amp30->GetXaxis()->SetTitle("x [mm]");
	g_UMass_amp30->GetXaxis()->SetTitle("x [mm]");


	gStyle->SetOptStat(0);

	new TCanvas();
	g_INFN_amp30->Draw("APZ");
	g_UMass_amp30->Draw("PZ");

	new TCanvas();
	g_INFN_blumlein->Draw("APZ");
	g_UMass_blumlein->Draw("PZ");

	TF1* f_quadratic = new TF1("f_quadratic","[0]+[1]*x*x");

	new TCanvas();
	g_INFN_amp30->SetMarkerColor(kGreen+2);
	g_UMass_amp30->SetMarkerColor(kOrange);
	g_INFN_blumlein->SetMarkerColor(kBlue);
	g_UMass_blumlein->SetMarkerColor(kRed);
	g_INFN_amp30->Draw("APZ");
	g_UMass_amp30->Draw("PZ");
	g_INFN_blumlein->Draw("PZ");
	g_UMass_blumlein->Draw("PZ");
	g_INFN_blumlein->Fit(f_quadratic);
	g_INFN_amp30->Fit(f_quadratic);
	gPad->BuildLegend();


	h1_K1_Rm6p6->SetLineColor(kBlue);
	h1_K1_R0->SetLineColor(kBlack);
	h1_K1_R3p2->SetLineColor(kRed-7);
	h1_K1_R6p6->SetLineColor(kRed);
	h1_K3_R0->SetLineColor(kGreen);

	h1_K1_Rm6p6_ra->SetLineColor(kBlue);
	h1_K1_R0_ra->SetLineColor(kBlack);
	h1_K1_R3p2_ra->SetLineColor(kRed-7);
	h1_K1_R6p6_ra->SetLineColor(kRed);
	h1_K3_R0_ra->SetLineColor(kGreen);

	h1_K1_Rm6p6->SetLineWidth(2);
	h1_K1_R0->SetLineWidth(2);
	h1_K1_R3p2->SetLineWidth(2);
	h1_K1_R6p6->SetLineWidth(2);
	h1_K3_R0->SetLineWidth(2);

	h1_K1_Rm6p6_ra->SetLineWidth(2);
	h1_K1_R0_ra->SetLineWidth(2);
	h1_K1_R3p2_ra->SetLineWidth(2);
	h1_K1_R6p6_ra->SetLineWidth(2);
	h1_K3_R0_ra->SetLineWidth(2);
	
	double ylow = -50;
	double yhigh = 20;
	
	h1_K1_Rm6p6->GetXaxis()->SetRangeUser(-0.1,0.3);
	h1_K1_Rm6p6_ra->GetXaxis()->SetRangeUser(-0.1,0.3);
	h1_K1_Rm6p6->GetYaxis()->SetRangeUser(ylow,yhigh);
	h1_K1_Rm6p6_ra->GetYaxis()->SetRangeUser(ylow,yhigh);
	h1_K1_R0->GetXaxis()->SetRangeUser(-0.1,0.3);
	h1_K1_R0_ra->GetXaxis()->SetRangeUser(-0.1,0.3);
	h1_K1_R0->GetYaxis()->SetRangeUser(ylow,yhigh);
	h1_K1_R0_ra->GetYaxis()->SetRangeUser(ylow,yhigh);

	TLine* l30 = new TLine(0.03,ylow,0.03,yhigh);
	l30->SetLineWidth(2);



	h1_K1_Rm6p6_ra_norm->SetLineColor(kBlue);
	h1_K1_R0_ra_norm->SetLineColor(kBlack);
	h1_K1_R3p2_ra_norm->SetLineColor(kRed-7);
	h1_K1_R6p6_ra_norm->SetLineColor(kRed);
	h1_K3_R0_ra_norm->SetLineColor(kBlack);
	h1_kick1_R0_ra_norm->SetLineColor(kBlack);
	h1_kick1_R1_ra_norm->SetLineColor(kRed+2);

	h1_K1_Rm6p6_ra_norm->SetLineWidth(2);
	h1_K1_R0_ra_norm->SetLineWidth(2);
	h1_K1_R3p2_ra_norm->SetLineWidth(2);
	h1_K1_R6p6_ra_norm->SetLineWidth(2);
	h1_K3_R0_ra_norm->SetLineWidth(2);
	h1_kick1_R0_ra_norm->SetLineWidth(2);
	h1_kick1_R1_ra_norm->SetLineWidth(2);

	new TCanvas();
	h1_K1_Rm6p6_ra_norm->Draw("HIST");
	h1_K1_R0_ra_norm->Draw("HIST SAME");
	h1_K1_R3p2_ra_norm->Draw("HIST SAME");
	h1_K1_R6p6_ra_norm->Draw("HIST SAME");
	h1_kick1_R0_ra_norm->Draw("HIST SAME");
	h1_kick1_R1_ra_norm->Draw("HIST SAME");
	gPad->BuildLegend();
	gPad->SetGridy();



	new TCanvas();
	h1_K1_Rm6p6_ra->Draw("HIST");
	h1_K1_R0_ra->Draw("HIST SAME");
	h1_K1_R3p2_ra->Draw("HIST SAME");
	h1_K1_R6p6_ra->Draw("HIST SAME");
	gPad->BuildLegend();
	l30->Draw("SAME");
	gPad->SetGridx();
	gPad->SetGridy();
	
	new TCanvas();
	h1_K1_Rm6p6->Draw("HIST");
	h1_K1_R0->Draw("HIST SAME");
	h1_K1_R3p2->Draw("HIST SAME");
	h1_K1_R6p6->Draw("HIST SAME");
	gPad->BuildLegend();
	l30->Draw("SAME");
	gPad->SetGridx();
	gPad->SetGridy();

	new TCanvas();
	h1_K1_R0->Draw("HIST");
	h1_K1_R0_ra->Draw("HIST SAME");
	gPad->BuildLegend();
	l30->Draw("SAME");
	gPad->SetGridx();
	gPad->SetGridy();


	h1_kick1_R0_ra->Scale(126./129.);

	h1_kick1_R0_ra->GetXaxis()->SetTitle("Time [ms]");
	h1_K3_R0_ra->GetXaxis()->SetTitle("Time [ms]");
	h1_kick1_R0_ra->GetYaxis()->SetTitle("B field [mG]");
	h1_K3_R0_ra->GetYaxis()->SetTitle("B field [mG]");

	h1_kick1_R0_ra->SetLineColor(kBlue);
	h1_K3_R0_ra->SetLineColor(kRed);
	h1_kick1_R0_ra->SetLineWidth(2);
	h1_K3_R0_ra->SetLineWidth(2);

	h1_kick1_R0_ra->GetXaxis()->SetRangeUser(-0.6,2.0);
	h1_K3_R0_ra->GetXaxis()->SetRangeUser(-0.6,2.0);
	h1_kick1_R0_ra->GetYaxis()->SetRangeUser(-60,140);
	h1_K3_R0_ra->GetYaxis()->SetRangeUser(-60,140);

	TH1D* h1_kick1_R0_ra_zoom = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_zoom");
	TH1D* h1_K3_R0_ra_zoom = (TH1D*)h1_K3_R0_ra->Clone("h1_K3_R0_ra_zoom");

	TCanvas* can = new TCanvas("can","",1200,900);
	h1_kick1_R0_ra->Draw("HIST");
	h1_K3_R0_ra->Draw("HIST SAME");
	gPad->SetGridy();
	can->SetHighLightColor(kBlack);
	TPad* pad_zoom = new TPad("pad_zoom","",0.4,0.40,0.85,0.88);
	pad_zoom->SetBorderMode(1);
	pad_zoom->Draw();
	pad_zoom->cd();
	h1_kick1_R0_ra_zoom->GetXaxis()->SetRangeUser(0,0.7);
	h1_K3_R0_ra_zoom->GetXaxis()->SetRangeUser(0,0.7);
	h1_kick1_R0_ra_zoom->GetYaxis()->SetRangeUser(-30,20);
	h1_K3_R0_ra_zoom->GetYaxis()->SetRangeUser(-30,20);
	h1_kick1_R0_ra_zoom->Draw("HIST");
	h1_K3_R0_ra_zoom->Draw("HIST SAME");
	TLegend* leg = new TLegend(0.4,0.2,0.6,0.4);
	leg->AddEntry(h1_kick1_R0_ra_zoom,"INFN","L");
	leg->AddEntry(h1_K3_R0_ra_zoom,"UMass","L");
	leg->Draw();
	TLine* l30zoom = new TLine(0.03,-30,0.03,20);
	l30zoom->SetLineWidth(2);
	l30zoom->SetLineStyle(kDashed);
	l30zoom->Draw("SAME");



	new TCanvas();
	h1_kick1_R0_ra->Draw("HIST");
	h1_K1_R0->Scale(55./53.);
	h1_K1_R0->Draw("HIST SAME");
	h1_K3_R0->Draw("HIST SAME");
}
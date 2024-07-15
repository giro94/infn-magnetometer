#include "../analysis_tools.C"

void plot_INFN_UMass(){

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;

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

	TGraph* g_INFN_blumlein = new TGraph();
	TGraph* g_INFN_amp30 = new TGraph();
	TGraph* g_UMass_blumlein = new TGraph();
	TGraph* g_UMass_amp30 = new TGraph();

	TH1D* h1_K1_R0_ra = smoothing(h1_K1_R0,"");
	TH1D* h1_K1_R3p2_ra = smoothing(h1_K1_R3p2,"");
	TH1D* h1_K1_R6p6_ra = smoothing(h1_K1_R6p6,"");
	TH1D* h1_K1_Rm6p6_ra = smoothing(h1_K1_Rm6p6,"");
	TH1D* h1_K3_R0_ra = smoothing(h1_K3_R0,"");



	g_INFN_blumlein->AddPoint(0,h1_kick1_R0_ra->Interpolate(-0.32)/h1_kick1_R0_ra->Interpolate(-0.32));
	g_INFN_blumlein->AddPoint(17.5,h1_kick1_R1_ra->Interpolate(-0.32)/h1_kick1_R0_ra->Interpolate(-0.32));

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

	g_INFN_blumlein->GetXaxis()->SetLimits(-20,20);
	g_UMass_blumlein->GetXaxis()->SetLimits(-20,20);
	g_INFN_amp30->GetXaxis()->SetLimits(-20,20);
	g_UMass_amp30->GetXaxis()->SetLimits(-20,20);
	g_INFN_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_UMass_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_INFN_amp30->GetXaxis()->SetTitle("x [mm]");
	g_UMass_amp30->GetXaxis()->SetTitle("x [mm]");


	new TCanvas();
	g_INFN_amp30->Draw("AP");
	g_UMass_amp30->Draw("P");

	new TCanvas();
	g_INFN_blumlein->Draw("AP");
	g_UMass_blumlein->Draw("P");


	new TCanvas();
	g_INFN_amp30->Draw("AP");
	g_UMass_amp30->Draw("P");
	g_INFN_blumlein->Draw("P");
	g_UMass_blumlein->Draw("P");


	new TCanvas();
	TH1D* h1_K1_R0_ra->Draw("HIST");
	TH1D* h1_K1_R3p2_ra->Draw("HIST SAME");
	TH1D* h1_K1_R6p6_ra->Draw("HIST SAME");
	TH1D* h1_K1_Rm6p6_ra->Draw("HIST SAME");
	TH1D* h1_K3_R0_ra->Draw("HIST SAME");
}
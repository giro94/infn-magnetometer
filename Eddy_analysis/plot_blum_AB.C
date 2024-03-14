#include "../analysis_tools.C"


void plot_blum_AB(){

	TFile* f1 = TFile::Open("analysis/analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root");


	TGraph* g_trend_A = (TGraph*)f1->Get("g_trend_A");
	TGraph* g_trend_B = (TGraph*)f1->Get("g_trend_B");
	TGraph* g_trend_ABdiff = (TGraph*)f1->Get("g_trend_ABdiff");
	TGraphErrors* g_trend_ABsum = (TGraphErrors*)f1->Get("g_trend_ABsum");
	TGraphErrors* g_trend_blumlein = (TGraphErrors*)f1->Get("g_trend_blumlein");
	TGraphErrors* g_trend_blumAB = (TGraphErrors*)f1->Get("g_trend_blumAB");


	TH1D* h1_blum = new TH1D("h1_blum","Blumlein distribution;Blumlein [mV]",100,0,100);
	TH1D* h1_blum_norm = new TH1D("h1_blum_norm","Normalized blumlein distribution;Blumlein [mV]",100,0,100);


	//Find maximum for ABsum
	double max = 0;
	int maxi = 0;
	for (int i=0; i<g_trend_ABsum->GetN(); i++){
		if (g_trend_ABsum->GetPointY(i) > max){
			max = g_trend_ABsum->GetPointY(i);
			maxi = i;
		}
	}

	//normalize g_trend_blumAB to first point
	double first_blum = g_trend_blumlein->GetPointY(maxi);
	double first_blumAB = g_trend_blumAB->GetPointY(maxi);
	for (int i=0; i<g_trend_blumAB->GetN(); i++){
		g_trend_blumAB->SetPointY(i,g_trend_blumAB->GetPointY(i)*first_blum/first_blumAB);


		h1_blum->Fill(g_trend_blumlein->GetPointY(i));
		h1_blum_norm->Fill(g_trend_blumAB->GetPointY(i));
	}

	g_trend_blumlein->SetMarkerStyle(20);
	g_trend_blumAB->SetMarkerStyle(20);
	g_trend_ABsum->SetMarkerStyle(20);

	g_trend_blumlein->SetMarkerColor(kBlue);
	g_trend_blumAB->SetMarkerColor(kRed);
	g_trend_ABsum->SetMarkerColor(kBlack);

	g_trend_ABsum->SetTitle("A+B trend");
	g_trend_blumlein->SetTitle("Blumlein trend");
	g_trend_blumAB->SetTitle("Blumlein trend normalized");


	gStyle->SetOptStat(0);


	h1_blum->SetLineColor(kBlue);
	h1_blum_norm->SetLineColor(kRed);
	h1_blum->SetLineWidth(2);
	h1_blum_norm->SetLineWidth(2);

	new TCanvas();
	g_trend_blumlein->Draw("APZ");
	g_trend_ABsum->Draw("PZ");
	gPad->BuildLegend();

	new TCanvas();
	g_trend_blumlein->Draw("APZ");
	g_trend_blumAB->Draw("PZ");
	gPad->BuildLegend();

	new TCanvas();
	h1_blum->Draw("HIST");
	h1_blum_norm->Draw("HIST SAME");
	gPad->BuildLegend();

}
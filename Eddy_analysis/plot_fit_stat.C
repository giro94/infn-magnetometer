void plot_fit_stat(){

	vector<TString> filenames = {
		"fits/fitted_noSNR_analysis_EC_jan22_B5173_H15.root",
		"fits/fitted_noSNR_analysis_EC_jan23_B5173_H25_Q00.root",
		"fits/fitted_noSNR_analysis_EC_jan26_B5173_H25_Q130.root",
		"fits/fitted_noSNR_analysis_EC_jan28_B5173_H25_Q00.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct10_H20_nofilter_Bfield.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct10_H22p5_nofilter_Bfield.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield_k777.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct21_H25_B100_k777.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct22_H25_B100_k777.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct24_H25_B100_k777.root",
		"fits/fitted_noSNR_analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root",
		"fits/fitted_noSNR_analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"fits/fitted_noSNR_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"fits/fitted_noSNR_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
	};

	vector<TString> filenames2 = {
		"analysis/analysis_EC_jan22_B5173_H15.root",
		"analysis/analysis_EC_jan23_B5173_H25_Q00.root",
		"analysis/analysis_EC_jan26_B5173_H25_Q130.root",
		"analysis/analysis_EC_jan28_B5173_H25_Q00.root",
		"analysis/analysis_SD_R0_eddy_oct10_H20_nofilter_Bfield.root",
		"analysis/analysis_SD_R0_eddy_oct10_H22p5_nofilter_Bfield.root",
		"analysis/analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield_k777.root",
		"analysis/analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield.root",
		"analysis/analysis_SD_R0_eddy_oct21_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct22_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct24_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root",
		"analysis/analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
	};

	new TCanvas();
	gStyle->SetPalette(kRainBow);

	TGraphErrors* g_exp_amp = new TGraphErrors();
	TGraphErrors* g_exp_tau = new TGraphErrors();
	for (int i=0; i<filenames.size(); i++){
		TFile* f_in = TFile::Open(filenames[i]);

		TH1D* h1_kick1_fit_exp = (TH1D*)f_in->Get("kick1_fit_exp17");
		h1_kick1_fit_exp->DrawCopy(i==0?"HIST PLC":"HIST SAME PLC");

		double a = h1_kick1_fit_exp->GetFunction("f_exp_17")->GetParameter(0);
		double da = h1_kick1_fit_exp->GetFunction("f_exp_17")->GetParError(0);
		double tau = h1_kick1_fit_exp->GetFunction("f_exp_17")->GetParameter(1);
		double dtau = h1_kick1_fit_exp->GetFunction("f_exp_17")->GetParError(1);

		g_exp_amp->SetPoint(i,i,a);
		g_exp_amp->SetPointError(i,0,da);
		g_exp_tau->SetPoint(i,i,1000*tau);
		g_exp_tau->SetPointError(i,0,1000*dtau);

		f_in->Close();
	}


	TGraph* g_AB_blum_all = new TGraph();
	TGraph** g_AB_blum = new TGraph*[filenames2.size()];
	TGraph** g_A = new TGraph*[filenames2.size()];
	TGraph** g_B = new TGraph*[filenames2.size()];
	TGraph** g_blum = new TGraph*[filenames2.size()];

	TGraph* g_blum_all = new TGraph();
	TGraph* g_AB_all = new TGraph();
	TGraph* g_ABratio_all = new TGraph();
	TH1D* h1_ABblumratio = new TH1D("h1_ABblumratio","",10,3,6);

	for (int i=0; i<filenames2.size(); i++){
		TFile* f_in = TFile::Open(filenames2[i]);

		g_AB_blum[i] = (TGraph*)f_in->Get("g_correlation_AB_blumlein");
		g_A[i] = (TGraph*)f_in->Get("g_trend_A");
		g_B[i] = (TGraph*)f_in->Get("g_trend_B");
		g_blum[i] = (TGraph*)f_in->Get("g_trend_blumlein");
		for (int j=0; j<g_AB_blum[i]->GetN(); j++){
			g_AB_blum_all->AddPoint(g_AB_blum[i]->GetPointX(j),g_AB_blum[i]->GetPointY(j));
		}
		f_in->Close();

		TFitResultPtr fitA = g_A[i]->Fit("pol0","S");
		TFitResultPtr fitB = g_B[i]->Fit("pol0","S");
		TFitResultPtr fitblum = g_blum[i]->Fit("pol0","S");
		double AB = fitA->Parameter(0)+fitB->Parameter(0);
		double blum = abs(fitblum->Parameter(0));
		double ABblumratio = blum/AB;

		g_blum_all->SetPoint(i,i,blum);
		g_AB_all->SetPoint(i,i,AB);
		g_ABratio_all->SetPoint(i,i,ABblumratio);
		h1_ABblumratio->Fill(ABblumratio);
	}

	new TCanvas();
	g_exp_amp->SetMarkerStyle(20);
	g_exp_amp->GetYaxis()->SetTitle("Amplitude [mG]");
	g_exp_amp->Draw("APL");
	//g_exp_amp->Fit("pol0");


	new TCanvas();
	g_exp_tau->SetMarkerStyle(20);
	g_exp_tau->GetYaxis()->SetTitle("Lifetime [#mus]");
	g_exp_tau->Draw("APL");
	//g_exp_tau->Fit("pol0");


	new TCanvas();
	g_AB_blum_all->SetMarkerStyle(20);
	g_AB_blum_all->Draw("AP");


	for (int i=0; i<filenames2.size(); i++){
		g_AB_blum[i]->SetMarkerStyle(20);
		g_AB_blum[i]->SetMarkerColor(1+(i%8));
		g_AB_blum[i]->Draw(i==0?"AP":"P");
	}


	new TCanvas();
	g_blum_all->SetMarkerStyle(20);
	g_blum_all->SetTitle("Blumlein");
	g_blum_all->Draw("APL");

	new TCanvas();
	g_AB_all->SetMarkerStyle(20);
	g_AB_all->SetTitle("A+B");
	g_AB_all->Draw("APL");

	new TCanvas();
	g_ABratio_all->SetMarkerStyle(20);
	g_ABratio_all->SetTitle("Blumlein / A+B");
	g_ABratio_all->Draw("APL");


	new TCanvas();
	h1_ABblumratio->SetTitle("Blumlein / A+B");
	h1_ABblumratio->Draw("HIST");

}
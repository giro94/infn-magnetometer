#include "../analysis_tools.C"

void plot_fit_stat(){



	vector<TString> filenames = {
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct21_H25_B100_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct22_H25_B100_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct24_H25_B100_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_dec15_H25_B100_k777.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_dec17_H25_B100_k777_long.root",
		"fits/fitted_calibrated_analysis_EC_jan22_B5173_H15.root",
		"fits/fitted_calibrated_analysis_EC_jan23_B5173_H25_Q00.root",
		"fits/fitted_calibrated_analysis_EC_jan26_B5173_H25_Q130.root",
		"fits/fitted_calibrated_analysis_EC_jan28_B5173_H25_Q00.root",
		"fits/fitted_calibrated_analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"fits/fitted_calibrated_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"fits/fitted_calibrated_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
	};

	vector<TString> filenames2 = {
		"analysis/analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield_k777.root",
		"analysis/analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield.root",
		"analysis/analysis_SD_R0_eddy_oct21_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct22_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct24_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_dec15_H25_B100_k777.root",
		"analysis/analysis_SD_R0_eddy_dec17_H25_B100_k777_long.root",
		"analysis/analysis_EC_jan22_B5173_H15.root",
		"analysis/analysis_EC_jan23_B5173_H25_Q00.root",
		"analysis/analysis_EC_jan26_B5173_H25_Q130.root",
		"analysis/analysis_EC_jan28_B5173_H25_Q00.root",
		"analysis/analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
	};





//noB
/*
	vector<TString> filenames = {
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct14_H25_nofilter_B0.root",
		"fits/fitted_calibrated_analysis_SD_R0_eddy_oct23_H25_B0_k777.root",
		"fits/fitted_calibrated_analysis_SD_R1_eddy_oct4_H0_nofilter.root",
	};

	vector<TString> filenames2 = {
		"analysis/analysis_SD_R0_eddy_oct14_H25_nofilter_B0.root",
		"analysis/analysis_SD_R0_eddy_oct23_H25_B0_k777.root",
		"analysis/analysis_SD_R1_eddy_oct4_H0_nofilter.root",
	};
*/
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

		g_exp_amp->SetPoint(i,i+1,a);
		g_exp_amp->SetPointError(i,0,da);
		g_exp_tau->SetPoint(i,i+1,1000*tau);
		g_exp_tau->SetPointError(i,0,1000*dtau);

		f_in->Close();
	}


	TGraphErrors* g_AB_blum_all = new TGraphErrors();
	TGraphErrors** g_AB = new TGraphErrors*[filenames2.size()];
	TGraphErrors** g_blum = new TGraphErrors*[filenames2.size()];
	TGraphErrors** g_blumAB = new TGraphErrors*[filenames2.size()];

	TGraphErrors* g_blum_all = new TGraphErrors();
	TGraphErrors* g_AB_all = new TGraphErrors();
	TGraphErrors* g_ABratio_all = new TGraphErrors();
	TH1D* h1_ABblumratio = new TH1D("h1_ABblumratio","",10,3,6);

	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");
	TGraphErrors* g_blumlein_ABnormalized = new TGraphErrors();
	TGraphErrors* g_blumlein_calibrated = new TGraphErrors();

	TH1D** h1_kick1_norm = new TH1D*[filenames2.size()];
	TH1D** h1_kick1_norm_ra = new TH1D*[filenames2.size()];

	for (int i=0; i<filenames2.size(); i++){
		TFile* f_in = TFile::Open(filenames2[i]);


		h1_kick1_norm[i] = ((TProfile*)f_in->Get("trace_kick1_ABnormalized"))->ProjectionX();
		cleanTrace(h1_kick1_norm[i],-50);
		h1_kick1_norm_ra[i] = runningAverage_5_10_15(h1_kick1_norm[i]);
		h1_kick1_norm_ra[i]->SetTitle(filenames2[i]);

		TH1D* kick1_calib = ((TProfile*)f_in->Get("trace_kick1_calibrated"))->ProjectionX();
		cleanTrace(kick1_calib,-150);
		TH1D* kick1_calib_ra = runningAverage_5_10_15(kick1_calib);

		cout<<h1_kick1_norm_ra[i]->Interpolate(0.03)<<"\n";

		new TCanvas();
		f_blumlein->SetParameters(50.0,-0.3,-1000.0);
		//Negative peak
		if (h1_kick1_norm_ra[i]->Interpolate(-0.3) < 0){
			f_blumlein->SetParameters(-50.0,-0.3,+1000.0);
		}
		TFitResultPtr fit_blumlein_norm = h1_kick1_norm_ra[i]->Fit("f_blumlein","QS+","",-0.4,-0.2);
		double fitted_blumlein_norm=0;
		double fitted_blumlein_norm_err=0;
		if (fit_blumlein_norm>=0){
			fitted_blumlein_norm = abs(fit_blumlein_norm->Parameter(0));
			fitted_blumlein_norm_err = fit_blumlein_norm->ParError(0);
		}
		g_blumlein_ABnormalized->SetPoint(i,i+1,fitted_blumlein_norm);
		g_blumlein_ABnormalized->SetPointError(i,0,fitted_blumlein_norm_err);

		f_blumlein->SetParameters(140.0,-0.3,-1000.0);
		TFitResultPtr fit_blumlein_calib = kick1_calib_ra->Fit("f_blumlein","N0QS+","",-0.4,-0.2);
		double fitted_blumlein_calib=0;
		double fitted_blumlein_calib_err=0;
		if (fit_blumlein_calib>=0){
			fitted_blumlein_calib = fit_blumlein_calib->Parameter(0);
			fitted_blumlein_calib_err = fit_blumlein_calib->ParError(0);
		}
		g_blumlein_calibrated->SetPoint(i,i+1,fitted_blumlein_calib);
		g_blumlein_calibrated->SetPointError(i,0,fitted_blumlein_calib_err);



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
	g_blumlein_ABnormalized->SetMarkerStyle(20);
	g_blumlein_ABnormalized->SetTitle("ABnormalized blumlein");
	g_blumlein_ABnormalized->Draw("APL");

	new TCanvas();
	g_blumlein_calibrated->SetMarkerStyle(20);
	g_blumlein_calibrated->SetTitle("Calibrated blumlein");
	g_blumlein_calibrated->Draw("APL");



	for (int i=0; i<filenames2.size(); i++){
		cout<<"file :"<<filenames2[i]<<"\n";
	}


	cout<<"Blumleins:\n";
	for (int i=0; i<filenames2.size(); i++){
		cout<<g_blumlein_ABnormalized->GetPointY(i)<<"\t"<<g_blumlein_ABnormalized->GetErrorY(i)<<"\n";
	}




}
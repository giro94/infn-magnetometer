#include "../analysis_tools.C"

void plot_calibrations(){


	vector<TString> filenames = {
		"analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield_k777.root",
		"analysis_SD_R0_eddy_oct17_H25_nofilter_Bfield.root",
		"analysis_SD_R0_eddy_oct21_H25_B100_k777.root",
		"analysis_SD_R0_eddy_oct22_H25_B100_k777.root",
		"analysis_SD_R0_eddy_oct24_H25_B100_k777.root",
		"analysis_SD_R0_eddy_oct27-29_H25_B100_k777.root",
		"analysis_SD_R0_eddy_dec15_H25_B100_k777.root",
		"analysis_SD_R0_eddy_dec17_H25_B100_k777_long.root",
		"analysis_EC_jan22_B5173_H15.root",
		"analysis_EC_jan23_B5173_H25_Q00.root",
		"analysis_EC_jan26_B5173_H25_Q130.root",
		"analysis_EC_jan28_B5173_H25_Q00.root",
		"analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
	};

	vector<TString> labels = {
		"Oct 17",
		"Oct 18",
		"Oct 21",
		"Oct 22",
		"Oct 24",
		"Oct 27",
		"Dec 15",
		"Dec 17",
		"Jan 22",
		"Jan 23",
		"Jan 26",
		"Jan 28",
		"Oct 5 (R1)",
		"Oct 8 (R1)",
		"Oct 9 (R1)",
	};


	vector<double> blumlein_raw = {
		52.746,
		51.998,
		52.116,
		52.581,
		51.611,
		52.610,
		58.128,
		58.720,
		47.206,
		43.378,
		62.386,
		43.617,
		61.243,
		62.567,
		62.743,
	};

	vector<double> blumlein_raw_error = {
		0.037,
		0.048,
		0.050,
		0.041,
		0.075,
		0.021,
		0.023,
		0.045,
		0.075,
		0.029,
		0.050,
		0.019,
		0.035,
		0.048,
		0.029,
	};

	vector<double> transient_30us = {
		6.6899,
		6.05684,
		6.58934,
		6.75209,
		6.49293,
		6.43074,
		7.93197,
		7.43393,
		6.22042,
		5.83125,
		7.88214,
		5.69761,
		13.5763,
		14.124,
		14.197,
	};

	vector<double> transient_30us_error = {
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
	};

	vector<TString> filenames_ramp = {
		"output_Rampdown_R0_H25_oct19_B100%-25%.root",
		"output_Rampdown_R0_H25_oct19_B100%-25%.root",
		"output_Rampup_R0_H25_oct24_B0-100.root",
		"output_Rampup_R0_H25_oct24_B0-100.root",
		"output_Rampup_R0_H25_oct24_B0-100.root",
		"output_Rampup_R0_H25_oct24_B0-100.root",
		"output_RampDown_R0_H25_dec18.root",
		"output_RampDown_R0_H25_dec18.root",
		"output_Ramp_jan18_0to5175.root",
		"output_Ramp_jan29_H25Q00_5173to0.root",
		"output_Ramp_jan26_H25Q130_5173to2000to5173.root",
		"output_Ramp_jan29_H25Q00_5173to0.root",
		"output_Rampup_R1_H5_oct5.root",
		"output_Rampup_R1_H5_oct5.root",
		"output_Rampup_R1_H5_oct5.root",
	};

	vector<int> ramp_idx = {
		0,
		0,
		1,
		1,
		1,
		1,
		2,
		2,
		3,
		4,
		5,
		4,
		6,
		6,
		6
	};

	vector<double> ramp_calibration = {
		0.429,
		0.429,
		0.406,
		0.406,
		0.406,
		0.406,
		0.472,
		0.472,
		0.342,
		0.348,
		0.489,
		0.348,
		0.442,
		0.442,
		0.442,
	};

	vector<double> ramp_calibration_error = {
		0.019,
		0.019,
		0.013,
		0.013,
		0.013,
		0.013,
		0.015,
		0.015,
		0.012,
		0.008,
		0.010,
		0.008,
		0.012,
		0.012,
		0.012,
	};

	vector<double> ramp_calibration_HWPfactor = {
		1.000,
		1.000,
		1.000,
		1.000,
		1.000,
		1.000,
		1.000,
		1.000,
		1.137,
		1.000,
		1.000,
		1.000,
		0.889,
		0.889,
		0.889,
	};






	vector<double> scaled_ramp_calibration;
	vector<double> scaled_ramp_calibration_error;
	for (int i=0; i<ramp_calibration_HWPfactor.size(); i++){
		scaled_ramp_calibration.push_back(ramp_calibration[i]*ramp_calibration_HWPfactor[i]);
		scaled_ramp_calibration_error.push_back(ramp_calibration_error[i]*ramp_calibration_HWPfactor[i]);
	}

	vector<double> calibrated_blumlein;
	vector<double> calibrated_blumlein_error;
	for (int i=0; i<blumlein_raw.size(); i++){
		calibrated_blumlein.push_back(blumlein_raw[i]/scaled_ramp_calibration[i]);
		double ramp_err = scaled_ramp_calibration_error[i]/scaled_ramp_calibration[i];
		double blum_err = blumlein_raw_error[i]/blumlein_raw[i];
		calibrated_blumlein_error.push_back(calibrated_blumlein[i]*sqrt(ramp_err*ramp_err + blum_err*blum_err));
	}

	vector<double> calibrated_transient;
	vector<double> calibrated_transient_error;
	for (int i=0; i<transient_30us.size(); i++){
		calibrated_transient.push_back(transient_30us[i]/scaled_ramp_calibration[i]);
		double ramp_err = scaled_ramp_calibration_error[i]/scaled_ramp_calibration[i];
		double trans_err = transient_30us_error[i]/transient_30us[i];
		calibrated_transient_error.push_back(calibrated_transient[i]*sqrt(ramp_err*ramp_err + trans_err*trans_err));
	}

	vector<double> blumlein_transient_ratio;
	for (int i=0; i<calibrated_transient.size(); i++){
		blumlein_transient_ratio.push_back(calibrated_transient[i]/calibrated_blumlein[i]);
	}





	//hists

	int Npoints = calibrated_blumlein.size();

	TH1D** h1_trace_blumlein = new TH1D* [Npoints];
	TH1D** h1_trace_transient = new TH1D* [Npoints];
	TH1D** h1_trace_blumlein_ra = new TH1D* [Npoints];
	TH1D** h1_trace_transient_ra = new TH1D* [Npoints];
	TH1D** h1_trace_uncalibrated = new TH1D* [Npoints];
	TH1D** h1_trace_uncalibrated_ra = new TH1D* [Npoints];
	TH1D** h1_trace_calibrated = new TH1D* [Npoints];
	TH1D** h1_trace_calibrated_ra = new TH1D* [Npoints];
	for (int i=0; i<Npoints; i++){
		cout<<filenames[i]<<"\n";
		TFile* fin = TFile::Open("analysis/"+filenames[i]);

		h1_trace_uncalibrated[i] = ((TProfile*)fin->Get("trace_kick1_ABnormalized"))->ProjectionX();
		if (h1_trace_uncalibrated[i]->Interpolate(-0.3) < 0){
			h1_trace_uncalibrated[i]->Scale(-1);
		}
		cleanTrace(h1_trace_uncalibrated[i],-200);
		h1_trace_uncalibrated_ra[i] = runningAverage_5_10_15(h1_trace_uncalibrated[i]);
		h1_trace_uncalibrated_ra[i]->SetTitle(filenames[i]);


		h1_trace_calibrated[i] = ((TProfile*)fin->Get("trace_kick1_ABnormalized"))->ProjectionX();
		if (h1_trace_calibrated[i]->Interpolate(-0.3) < 0){
			h1_trace_calibrated[i]->Scale(-1);
		}
		cleanTrace(h1_trace_calibrated[i],-200);
		h1_trace_calibrated_ra[i] = runningAverage_5_10_15(h1_trace_calibrated[i]);
		h1_trace_calibrated_ra[i]->SetTitle(filenames[i]);
		h1_trace_calibrated_ra[i]->Scale(1./scaled_ramp_calibration[i]);


		h1_trace_blumlein[i] = ((TProfile*)fin->Get("trace_kick1_ABnormalized"))->ProjectionX();
		if (h1_trace_blumlein[i]->Interpolate(-0.3) < 0){
			h1_trace_blumlein[i]->Scale(-1);
		}
		cleanTrace(h1_trace_blumlein[i],-200);
		h1_trace_blumlein_ra[i] = runningAverage_5_10_15(h1_trace_blumlein[i]);
		h1_trace_blumlein_ra[i]->SetTitle(filenames[i]);
		h1_trace_blumlein_ra[i]->Scale(1./scaled_ramp_calibration[i]);


		h1_trace_transient[i] = ((TProfile*)fin->Get("trace_kick1_ABnormalized"))->ProjectionX();
		if (h1_trace_transient[i]->Interpolate(-0.3) < 0){
			h1_trace_transient[i]->Scale(-1);
		}
		cleanTrace(h1_trace_transient[i],-200);
		h1_trace_transient_ra[i] = runningAverage_5_10_15(h1_trace_transient[i]);
		h1_trace_transient_ra[i]->SetTitle(filenames[i]);
		h1_trace_transient_ra[i]->Scale(1./scaled_ramp_calibration[i]);

	}


	TH1D* h1_calibrated_blumlein = new TH1D("h1_calibrated_blumlein","Calibrated blumlein;;Blumlein [mG]",Npoints,0,Npoints);
	TH1D* h1_calibrated_transient = new TH1D("h1_calibrated_transient","Calibrated transient (30 #mus);;Transient (30 #mus) [mG]",Npoints,0,Npoints);
	TH1D** h1_calibrated_blumlein_calibs = new TH1D* [7];
	TH1D** h1_calibrated_transient_calibs = new TH1D* [7];
	for (int i=0; i<7; i++){
		h1_calibrated_blumlein_calibs[i] = new TH1D(Form("h1_calibrated_blumlein_calibs%d",i),"Calibrated blumlein;;Blumlein [mG]",Npoints,0,Npoints);
		h1_calibrated_transient_calibs[i] = new TH1D(Form("h1_calibrated_transient_calibs%d",i),"Calibrated transient (30 #mus);;Transient (30 #mus) [mG]",Npoints,0,Npoints);
	}

	for (int i=0; i<Npoints; i++){
		cout<<"("<<labels[i]<<") Blumlein: "<<calibrated_blumlein[i]<<" +- "<<calibrated_blumlein_error[i]<<" Transient: "<<calibrated_transient[i]<<" +- "<<calibrated_transient_error[i]<<"\n";
		h1_calibrated_blumlein->SetBinContent(i+1,calibrated_blumlein[i]);
		h1_calibrated_blumlein->SetBinError(i+1,calibrated_blumlein_error[i]);
		h1_calibrated_blumlein->GetXaxis()->SetBinLabel(i+1,labels[i]);
		h1_calibrated_transient->SetBinContent(i+1,calibrated_transient[i]);
		h1_calibrated_transient->SetBinError(i+1,calibrated_transient_error[i]);
		h1_calibrated_transient->GetXaxis()->SetBinLabel(i+1,labels[i]);

		int calib_idx = ramp_idx[i];
		h1_calibrated_blumlein_calibs[calib_idx]->SetBinContent(i+1,calibrated_blumlein[i]);
		h1_calibrated_blumlein_calibs[calib_idx]->SetBinError(i+1,calibrated_blumlein_error[i]);
		h1_calibrated_blumlein_calibs[calib_idx]->GetXaxis()->SetBinLabel(i+1,labels[i]);
		h1_calibrated_transient_calibs[calib_idx]->SetBinContent(i+1,calibrated_transient[i]);
		h1_calibrated_transient_calibs[calib_idx]->SetBinError(i+1,calibrated_transient_error[i]);
		h1_calibrated_transient_calibs[calib_idx]->GetXaxis()->SetBinLabel(i+1,labels[i]);
	}

	
	TH1D* h1_blumlein_R0 = new TH1D("h1_blumlein_R0","Blumlein R0;Blumlein [mG]",35,100,170);
	TH1D* h1_blumlein_R1 = new TH1D("h1_blumlein_R1","Blumlein R1;Blumlein [mG]",35,100,170);
	TH1D* h1_transient_R0 = new TH1D("h1_transient_R0","Transient R0;Transient [mG]",50,0,50);
	TH1D* h1_transient_R1 = new TH1D("h1_transient_R1","Transient R1;Transient [mG]",50,0,50);
	TH1D* h1_ratio_R0 = new TH1D("h1_ratio_R0","Transient/Blumlein ratio R0;Transient/blumlein",25,0,0.3);
	TH1D* h1_ratio_R1 = new TH1D("h1_ratio_R1","Transient/Blumlein ratio R1;Transient/blumlein",25,0,0.3);
	
	for (int i=0; i<Npoints; i++){
		if (i<Npoints-3) {
			h1_blumlein_R0->Fill(calibrated_blumlein[i]);
			h1_transient_R0->Fill(calibrated_transient[i]);
			h1_ratio_R0->Fill(calibrated_transient[i]/calibrated_blumlein[i]);
		}
		else {
			h1_blumlein_R1->Fill(calibrated_blumlein[i]);
			h1_transient_R1->Fill(calibrated_transient[i]);
			h1_ratio_R1->Fill(calibrated_transient[i]/calibrated_blumlein[i]);
		}
	}




	//drawing

	gStyle->SetOptStat(0);

	vector<int> colors = {1,2,3,4,6,kOrange,kGreen+2};


	for (int i=0; i<7; i++){
		h1_calibrated_blumlein_calibs[i]->LabelsOption("v");
		h1_calibrated_blumlein_calibs[i]->SetMarkerStyle(20);
		h1_calibrated_blumlein_calibs[i]->SetMarkerSize(1.5);
		h1_calibrated_blumlein_calibs[i]->SetMarkerColor(colors[i]);
		h1_calibrated_blumlein_calibs[i]->SetLineColor(colors[i]);
		h1_calibrated_blumlein_calibs[i]->SetLineWidth(2);
		h1_calibrated_transient_calibs[i]->LabelsOption("v");
		h1_calibrated_transient_calibs[i]->SetMarkerStyle(20);
		h1_calibrated_transient_calibs[i]->SetMarkerSize(1.5);
		h1_calibrated_transient_calibs[i]->SetMarkerColor(colors[i]);
		h1_calibrated_transient_calibs[i]->SetLineColor(colors[i]);
		h1_calibrated_transient_calibs[i]->SetLineWidth(2);
	}

	new TCanvas();
	h1_calibrated_blumlein->LabelsOption("v");
	h1_calibrated_blumlein->SetMarkerStyle(20);
	h1_calibrated_blumlein->Draw("PEX0");
	for (int i=0; i<7; i++){
		h1_calibrated_blumlein_calibs[i]->Draw("PEX0 SAME");
	}
	TFitResultPtr fit_blum_R0 = h1_calibrated_blumlein->Fit("pol0","S+","",0.4,11.6);
	TFitResultPtr fit_blum_R1 = h1_calibrated_blumlein->Fit("pol0","S+","",12.4,14.6);

	new TCanvas();
	h1_calibrated_transient->LabelsOption("v");
	h1_calibrated_transient->SetMarkerStyle(20);
	h1_calibrated_transient->Draw("PEX0");
	for (int i=0; i<7; i++){
		h1_calibrated_transient_calibs[i]->Draw("PEX0 SAME");
	}
	TFitResultPtr fit_trans_R0 = h1_calibrated_transient->Fit("pol0","S+","",0.4,11.6);
	TFitResultPtr fit_trans_R1 = h1_calibrated_transient->Fit("pol0","S+","",12.4,14.6);


	new TCanvas("","",1400,1000);
	for (int i=0; i<Npoints; i++){
		h1_trace_uncalibrated_ra[i]->GetXaxis()->SetRangeUser(-0.6,0.7);
		h1_trace_uncalibrated_ra[i]->GetYaxis()->SetRangeUser(-40,80);
		h1_trace_uncalibrated_ra[i]->Draw(i==0?"HIST":"HIST SAME");
	}

	new TCanvas("","",1400,1000);
	for (int i=0; i<Npoints; i++){
		h1_trace_calibrated_ra[i]->GetXaxis()->SetRangeUser(-0.6,0.7);
		h1_trace_calibrated_ra[i]->GetYaxis()->SetRangeUser(-40,80);
		h1_trace_calibrated_ra[i]->Draw(i==0?"HIST":"HIST SAME");
	}


	new TCanvas();
	for (int i=0; i<Npoints; i++){
		h1_trace_blumlein_ra[i]->SetLineWidth(2);
		h1_trace_blumlein_ra[i]->SetLineColor(colors[ramp_idx[i]]);
		h1_trace_blumlein_ra[i]->GetXaxis()->SetRangeUser(-0.45,-0.2);
		h1_trace_blumlein_ra[i]->GetYaxis()->SetRangeUser(80,180);
		h1_trace_blumlein_ra[i]->Draw(i==0?"HIST":"HIST SAME");
	}

	new TCanvas();
	for (int i=0; i<Npoints; i++){
		h1_trace_transient_ra[i]->SetLineWidth(2);
		h1_trace_transient_ra[i]->SetLineColor(colors[ramp_idx[i]]);
		h1_trace_transient_ra[i]->GetXaxis()->SetRangeUser(0,0.3);
		h1_trace_transient_ra[i]->GetYaxis()->SetRangeUser(-80,10);
		h1_trace_transient_ra[i]->Draw(i==0?"HIST":"HIST SAME");
	}

	new TCanvas();
	for (int i=0; i<Npoints; i++){
		h1_trace_transient[i]->SetLineWidth(2);
		h1_trace_transient[i]->SetLineColor(colors[ramp_idx[i]]);
		h1_trace_transient[i]->GetXaxis()->SetRangeUser(0,0.3);
		h1_trace_transient[i]->GetYaxis()->SetRangeUser(-80,10);
		h1_trace_transient[i]->Draw(i==0?"HIST":"HIST SAME");
	}

	new TCanvas();
	h1_blumlein_R0->SetLineColor(kBlue);
	h1_blumlein_R1->SetLineColor(kRed);
	h1_blumlein_R0->SetLineWidth(2);
	h1_blumlein_R1->SetLineWidth(2);
	h1_blumlein_R0->Draw("HIST");
	h1_blumlein_R1->Draw("HIST SAME");

	new TCanvas();
	h1_transient_R0->SetLineColor(kBlue);
	h1_transient_R1->SetLineColor(kRed);
	h1_transient_R0->SetLineWidth(2);
	h1_transient_R1->SetLineWidth(2);
	h1_transient_R0->Draw("HIST");
	h1_transient_R1->Draw("HIST SAME");

	new TCanvas();
	h1_ratio_R0->SetLineColor(kBlue);
	h1_ratio_R1->SetLineColor(kRed);
	h1_ratio_R0->SetLineWidth(2);
	h1_ratio_R1->SetLineWidth(2);
	h1_ratio_R0->Draw("HIST");
	h1_ratio_R1->Draw("HIST SAME");


	double blum_R0_mean = fit_blum_R0->Parameter(0);
	double blum_R1_mean = fit_blum_R1->Parameter(0);
	double blum_R0_err = h1_blumlein_R0->GetRMS();
	double blum_R1_err = blum_R1_mean*(blum_R0_err/blum_R0_mean);

	double trans_R0_mean = fit_trans_R0->Parameter(0);
	double trans_R1_mean = fit_trans_R1->Parameter(0);
	double trans_R0_err = h1_transient_R0->GetRMS();
	double trans_R1_err = trans_R1_mean*(trans_R0_err/trans_R0_mean);

	cout<<"Blumlein\n";
	cout<<"R0: "<<blum_R0_mean<<" +- "<<blum_R0_err<<"\n";
	cout<<"R0: "<<blum_R1_mean<<" +- "<<blum_R1_err<<"\n";

	cout<<"Transient\n";
	cout<<"R0: "<<trans_R0_mean<<" +- "<<trans_R0_err<<"\n";
	cout<<"R0: "<<trans_R1_mean<<" +- "<<trans_R1_err<<"\n";

	TGraphErrors* g_blumlein = new TGraphErrors();
	TGraphErrors* g_transient = new TGraphErrors();

	g_blumlein->SetPoint(0,0,blum_R0_mean);
	g_blumlein->SetPointError(0,2,blum_R0_err);
	g_blumlein->SetPoint(1,17.5,blum_R1_mean);
	g_blumlein->SetPointError(1,2,blum_R1_err);

	g_transient->SetPoint(0,0,trans_R0_mean);
	g_transient->SetPointError(0,2,trans_R0_err);
	g_transient->SetPoint(1,17.5,trans_R1_mean);
	g_transient->SetPointError(1,2,trans_R1_err);

	new TCanvas();
	g_blumlein->SetMarkerColor(kBlue);
	g_blumlein->SetMarkerStyle(20);
	g_transient->SetMarkerColor(kRed);
	g_transient->SetMarkerStyle(34);
	g_blumlein->GetXaxis()->SetLimits(-45,45);
	g_blumlein->GetYaxis()->SetRangeUser(0,160);
	g_blumlein->Draw("APZ");
	g_transient->Draw("PZ");

}
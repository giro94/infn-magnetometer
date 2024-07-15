#include "../analysis_tools.C"




void plot_blumlein_R1(){


	vector<TString> filenames = {
		//"analysis/analysis_FD_R1_blum_oct8_H0_Bfield.root",
		//"analysis/analysis_SD_R1_eddy_oct4_H0_nofilter.root",
		"analysis/analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
		"analysis/analysis_SD_R1_eddy_oct7_H0_nofilter_Bfield_k194.root",
		"analysis/analysis_SD_R1_eddy_oct7_H0_nofilter_Bfield_k389.root",
		"analysis/analysis_SD_R1_eddy_oct6_H0_nofilter_Bfield_k583.root"
	};

	vector<TString> nicknames = {
		//"Fast Diodes (Oct 8)",
		//"Field OFF (Oct 4)",
		"Field ON (Oct 5)",
		"Field ON (Oct 8)",
		"Field ON (Oct 8)",
		"Kick 1.94 kV (Oct 7)",
		"Kick 3.89 kV (Oct 7)",
		"Kick 5.83 kV (Oct 6)"
	};
	
	vector<TString> filenames_fit = {
		//"fits/fitted_normalized_analysis_FD_R1_blum_oct8_H0_Bfield.root",
		//"fits/fitted_normalized_analysis_SD_R1_eddy_oct4_H0_nofilter.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct5_H0_nofilter_Bfield.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield_afterrampup.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct7_H0_nofilter_Bfield_k194.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct7_H0_nofilter_Bfield_k389.root",
		"fits/fitted_normalized_analysis_SD_R1_eddy_oct6_H0_nofilter_Bfield_k583.root"
	};

	vector<double> kick_strength = {
		//7.77,
		//7.77,
		7.77,
		7.77,
		7.77,
		1.94,
		3.89,
		5.83,
	};
	int Nfiles = filenames.size();
	TFile** f = new TFile*[Nfiles];
	TFile** f_fit = new TFile*[Nfiles];


	double mV_to_mG = 1./(0.442 * 0.91);


	vector<double> blum_avg;
	vector<double> blum_avg_err;
	vector<double> fit_amp;
	vector<double> fit_amp_err;
	TGraphErrors** g_ABsum = new TGraphErrors* [Nfiles];
	TGraphErrors** g_blum = new TGraphErrors* [Nfiles];
	TGraphErrors** g_blum_norm = new TGraphErrors* [Nfiles];
	TH1D** h1_fit_exp = new TH1D* [Nfiles];
	TH1D** h1_blum_norm = new TH1D* [Nfiles];
	for (int i=0; i<Nfiles; i++){

		f[i] = TFile::Open(filenames[i]);
		g_ABsum[i] = (TGraphErrors*)f[i]->Get("g_trend_ABsum");
		g_blum[i] = (TGraphErrors*)f[i]->Get("g_trend_blumlein");
		g_blum_norm[i] = (TGraphErrors*)f[i]->Get("g_trend_blumAB");

		f_fit[i] = TFile::Open(filenames_fit[i]);
		h1_fit_exp[i] = (TH1D*)f_fit[i]->Get("kick1_fit_exp");

		TString hname = Form("h1_blum_norm_%d",i);
		TString htitle = Form("Blumlein distribution (%s);Blumlein [mV]",nicknames[i].Data());
		h1_blum_norm[i] = new TH1D(hname,htitle,80,0,80);


		//Clean outliers
		double avg = g_blum_norm[i]->GetMean(2);
		double std_dev = g_blum_norm[i]->GetRMS(2);
		if (std_dev > 3) std_dev = 2;
		cout<<"File "<<i<<" mean "<<avg<<" std "<<std_dev<<"\n";
		for (int j=0; j<g_blum_norm[i]->GetN(); j++){
			double yval = g_blum_norm[i]->GetPointY(j);
			g_blum_norm[i]->SetPointX(j,j);
			if (abs(yval - avg) > 3*std_dev){
				g_blum_norm[i]->RemovePoint(j);
				j--;
			} else {
				h1_blum_norm[i]->Fill(yval);
			}
		}

		for (int j=0; j<g_ABsum[i]->GetN(); j++){
			g_ABsum[i]->SetPointX(j,j);
		}
		for (int j=0; j<g_blum[i]->GetN(); j++){
			g_blum[i]->SetPointX(j,j);
		}

		g_ABsum[i]->GetXaxis()->SetTimeDisplay(0);
		g_ABsum[i]->GetXaxis()->SetTimeFormat("");

		g_blum[i]->GetXaxis()->SetTimeDisplay(0);
		g_blum[i]->GetXaxis()->SetTimeFormat("");

		g_blum_norm[i]->GetXaxis()->SetTimeDisplay(0);
		g_blum_norm[i]->GetXaxis()->SetTimeFormat("");

		double new_avg = g_blum_norm[i]->GetMean(2);
		double new_std_dev = g_blum_norm[i]->GetRMS(2);
		double new_avg_err = new_std_dev/sqrt(g_blum_norm[i]->GetN());
		blum_avg.push_back(new_avg);
		blum_avg_err.push_back(new_avg_err);

		g_ABsum[i]->SetMarkerStyle(20);
		g_blum[i]->SetMarkerStyle(20);
		g_blum_norm[i]->SetMarkerStyle(20);

		g_ABsum[i]->SetMarkerColor(i+1);
		g_blum[i]->SetMarkerColor(i+1);
		g_blum_norm[i]->SetMarkerColor(i+1);

		h1_fit_exp[i]->SetLineWidth(2);;
		h1_fit_exp[i]->SetLineColor(i+1);
		h1_fit_exp[i]->GetXaxis()->SetRangeUser(-0.6,0.6);
		h1_fit_exp[i]->GetYaxis()->SetRangeUser(-40,70);
		h1_fit_exp[i]->SetTitle(Form("Trace (12 V) (%s);Time [ms];Trace [mV]",nicknames[i].Data()));

		fit_amp.push_back(h1_fit_exp[i]->GetFunction("f_exp")->GetParameter(0));
		fit_amp_err.push_back(h1_fit_exp[i]->GetFunction("f_exp")->GetParError(0));

		h1_blum_norm[i]->SetLineWidth(2);;
		h1_blum_norm[i]->SetLineColor(i+1);

		g_ABsum[i]->SetTitle(Form("A+B (%s);Time;A+B [mV]",nicknames[i].Data()));
		g_blum[i]->SetTitle(Form("Blumlein (raw) (%s);Time;Blumlein [mV]",nicknames[i].Data()));
		g_blum_norm[i]->SetTitle(Form("Blumlein (12 V) (%s);Time;Blumlein [mV]",nicknames[i].Data()));

	}



	for (int i=0; i<Nfiles; i++){
		cout<<"File "<<filenames[i]<<"\n";
		cout<<"Nick "<<nicknames[i]<<"\n";
		cout<<"Blumlein: "<<blum_avg[i]<<" +- "<<blum_avg_err[i]<<" mV\n";
	}

	gStyle->SetOptStat(0);

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_blum_norm[i]->GetYaxis()->SetRangeUser(0,70);
		g_blum_norm[i]->Draw(i==0?"APZ":"PZ");
	}
	gPad->BuildLegend();

	TCanvas* can = new TCanvas();
	can->Divide(3,1);
	can->cd(1);
	for (int i=0; i<Nfiles; i++){
		g_ABsum[i]->GetYaxis()->SetRangeUser(0,15);
		g_ABsum[i]->Draw(i==0?"APZ":"PZ");
	}
	gPad->BuildLegend();
	can->cd(2);
	for (int i=0; i<Nfiles; i++){
		g_blum[i]->GetYaxis()->SetRangeUser(0,70);
		g_blum[i]->Draw(i==0?"APZ":"PZ");
	}
	gPad->BuildLegend();
	can->cd(3);
	for (int i=0; i<Nfiles; i++){
		g_blum_norm[i]->GetYaxis()->SetRangeUser(0,70);
		g_blum_norm[i]->Draw(i==0?"APZ":"PZ");
	}
	gPad->BuildLegend();

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		h1_blum_norm[i]->Draw(i==0?"HIST":"HIST SAME");
	}
	gPad->BuildLegend();


	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		h1_fit_exp[i]->Draw(i==0?"HIST":"HIST SAME");
	}
	gPad->BuildLegend();
	for (int i=0; i<Nfiles; i++){
		h1_fit_exp[i]->GetFunction("f_exp")->Draw("SAME");
	}

	TF1* f_pol1_zero_blum = new TF1("f_pol1_zero_blum","[0]*x");
	TF1* f_pol1_zero_amp = new TF1("f_pol1_zero_amp","[0]*x");

	TGraphErrors* g_kick_blumlein = new TGraphErrors(6);
	for (int i=0; i<Nfiles; i++){
		g_kick_blumlein->SetPoint(i,kick_strength[i],blum_avg[i]*mV_to_mG);
		g_kick_blumlein->SetPointError(i,0,blum_avg_err[i]*mV_to_mG);
	}
	g_kick_blumlein->SetMarkerStyle(20);
	g_kick_blumlein->GetXaxis()->SetLimits(0,9);
	g_kick_blumlein->GetYaxis()->SetRangeUser(0,200);
	g_kick_blumlein->GetXaxis()->SetTitle("Kicker strength");
	g_kick_blumlein->GetYaxis()->SetTitle("Blumlein [mG]");
	g_kick_blumlein->SetTitle("Blumlein amplitude [mG]");
	new TCanvas();
	g_kick_blumlein->Draw("APZ");
	g_kick_blumlein->Fit("pol1");


	TGraphErrors* g_kick_amplitude = new TGraphErrors(6);
	for (int i=0; i<Nfiles; i++){
		g_kick_amplitude->SetPoint(i,kick_strength[i],fit_amp[i]*mV_to_mG);
		g_kick_amplitude->SetPointError(i,0,fit_amp_err[i]*mV_to_mG);
	}
	g_kick_amplitude->SetMarkerStyle(22);
	g_kick_amplitude->SetMarkerColor(kRed);
	g_kick_amplitude->GetXaxis()->SetLimits(0,9);
	g_kick_amplitude->GetYaxis()->SetRangeUser(0,200);
	g_kick_amplitude->GetXaxis()->SetTitle("Kicker strength");
	g_kick_amplitude->GetYaxis()->SetTitle("EC Amplitude [mG]");
	g_kick_amplitude->SetTitle("Transient amplitude [mG]");
	new TCanvas();
	g_kick_amplitude->Draw("APZ");
	g_kick_amplitude->Fit("pol1");


	new TCanvas();
	g_kick_blumlein->Draw("APZ");
	g_kick_blumlein->Fit("f_pol1_zero_blum");
	g_kick_amplitude->Draw("PZ");
	g_kick_amplitude->Fit("f_pol1_zero_amp");
	gPad->BuildLegend();

}
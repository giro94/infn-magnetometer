#include "../analysis_tools.C"

void plot_Bvalues(){

	TFile* f0_1 = TFile::Open("analysis/analysis_SD_R0_eddy_oct14_H25_nofilter_B0.root");
	TFile* f0_2 = TFile::Open("analysis/analysis_SD_R0_eddy_oct23_H25_B0_k777.root");

	TFile* f1_p = TFile::Open("analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f1_n = TFile::Open("analysis/analysis_EC_jan28_B5173_H25_Q00.root");
	TFile* f1_0 = TFile::Open("analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");

	TFile* f2 = TFile::Open("analysis/analysis_EC_jan21_B3043_H24.root");

	TFile* f3 = TFile::Open("analysis/analysis_EC_jan19_B4353_H23.root");

	double blum_norm_x = -0.323;
	double blum_norm_y = 129.;
	double blum_norm_x_trace = 4.769;
	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1_kick1_B0_1 = ((TProfile*)f0_1->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B0_2 = ((TProfile*)f0_2->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B1_p = ((TProfile*)f1_p->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_B1_n = ((TProfile*)f1_n->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_B1_0 = ((TProfile*)f1_0->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_B2 = ((TProfile*)f2->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_B3 = ((TProfile*)f3->Get("trace_kick1_calibrated"))->ProjectionX();

	TH1D* h1_kick1_B1 = (TH1D*)h1_kick1_B1_p->Clone("h1_kick1_B1");
	h1_kick1_B1->Scale(0.5);
	h1_kick1_B1->Add(h1_kick1_B1_n,-0.5);
	h1_kick1_B1->Add(h1_kick1_B1_0,0.2);

	TH1D* h1_kick1_B1_ra = smoothing(h1_kick1_B1,"");
	double b1_norm = blum_norm_y/h1_kick1_B1_ra->Interpolate(blum_norm_x);
	h1_kick1_B1_ra->Scale(b1_norm);
	h1_kick1_B1->Scale(b1_norm);


	TH1D* h1_kick1_B0a = (TH1D*)h1_kick1_B0_1->Clone("h1_kick1_B0a");
	h1_kick1_B0a->Add(h1_kick1_B0_2);
	h1_kick1_B0a->Scale(0.5);



	TH1D* h1_kick1_B0_1_ra = smoothing(h1_kick1_B0_1,"");
	TH1D* h1_kick1_B0_2_ra = smoothing(h1_kick1_B0_2,"");
	TH1D* h1_kick1_B0a_ra = smoothing(h1_kick1_B0a,"");
	TH1D* h1_kick1_B1_p_ra = smoothing(h1_kick1_B1_p,"");
	TH1D* h1_kick1_B1_n_ra = smoothing(h1_kick1_B1_n,"");
	TH1D* h1_kick1_B1_0_ra = smoothing(h1_kick1_B1_0,"");
	TH1D* h1_kick1_B2_ra = smoothing(h1_kick1_B2,"");
	TH1D* h1_kick1_B3_ra = smoothing(h1_kick1_B3,"");

	TH1D* h1_kick1_B0 = (TH1D*)h1_kick1_B0_1->Clone("h1_kick1_B0");
	TH1D* h1_kick1_B0_ra = (TH1D*)h1_kick1_B0_1_ra->Clone("h1_kick1_B0_ra");


	//////////// Re-baseline

	TF1* f_baseline = new TF1("f_baseline","[0]",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B0_1_ra = h1_kick1_B0_1_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B0_2_ra = h1_kick1_B0_2_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B0a_ra = h1_kick1_B0a_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B1_p_ra = h1_kick1_B1_p_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B1_n_ra = h1_kick1_B1_n_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B1_0_ra = h1_kick1_B1_0_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B2_ra = h1_kick1_B2_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B3_ra = h1_kick1_B3_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	TFitResultPtr fit_h1_kick1_B0_ra = h1_kick1_B0_ra->Fit(f_baseline,"QNS","",-1,-0.6);
	
	addScalar(h1_kick1_B0_1_ra,-fit_h1_kick1_B0_1_ra->Parameter(0));
	addScalar(h1_kick1_B0_2_ra,-fit_h1_kick1_B0_2_ra->Parameter(0));
	addScalar(h1_kick1_B0a_ra,-fit_h1_kick1_B0a_ra->Parameter(0));
	addScalar(h1_kick1_B1_p_ra,-fit_h1_kick1_B1_p_ra->Parameter(0));
	addScalar(h1_kick1_B1_n_ra,-fit_h1_kick1_B1_n_ra->Parameter(0));
	addScalar(h1_kick1_B1_0_ra,-fit_h1_kick1_B1_0_ra->Parameter(0));
	addScalar(h1_kick1_B2_ra,-fit_h1_kick1_B2_ra->Parameter(0));
	addScalar(h1_kick1_B3_ra,-fit_h1_kick1_B3_ra->Parameter(0));
	addScalar(h1_kick1_B0_ra,-fit_h1_kick1_B0_ra->Parameter(0));

	h1_kick1_B0_1_ra->Scale(blum_norm_y/h1_kick1_B0_1_ra->Interpolate(blum_norm_x));
	h1_kick1_B0_2_ra->Scale(blum_norm_y/h1_kick1_B0_2_ra->Interpolate(blum_norm_x));
	h1_kick1_B0a_ra->Scale(blum_norm_y/h1_kick1_B0a_ra->Interpolate(blum_norm_x));
	h1_kick1_B1_p_ra->Scale(blum_norm_y/h1_kick1_B1_p_ra->Interpolate(blum_norm_x));
	h1_kick1_B1_n_ra->Scale(blum_norm_y/h1_kick1_B1_n_ra->Interpolate(blum_norm_x));
	h1_kick1_B1_0_ra->Scale(blum_norm_y/h1_kick1_B1_0_ra->Interpolate(blum_norm_x));
	h1_kick1_B2_ra->Scale(blum_norm_y/h1_kick1_B2_ra->Interpolate(blum_norm_x));
	h1_kick1_B3_ra->Scale(blum_norm_y/h1_kick1_B3_ra->Interpolate(blum_norm_x));
	h1_kick1_B0_ra->Scale(blum_norm_y/h1_kick1_B0_ra->Interpolate(blum_norm_x));


	//////////////////////////////////////

	double xmin = -1;
	double xmax = 8;
	int nbins = 8928;

	TH1D* h1_kick1_B0_ra_resampled = new TH1D("h1_kick1_B0_ra_resampled","",nbins,xmin,xmax);
	for (int bx=1; bx<=h1_kick1_B0_ra_resampled->GetNbinsX(); bx++){
		double bin_center = h1_kick1_B0_ra_resampled->GetBinCenter(bx);
		double yval = h1_kick1_B0_ra->Interpolate(bin_center);
		h1_kick1_B0_ra_resampled->SetBinContent(bx,yval);
	}

	TH1D* h1_kick1_B1_B0_ra = (TH1D*)h1_kick1_B1_ra->Clone("h1_kick1_B1_B0_ra");
	h1_kick1_B1_B0_ra->Add(h1_kick1_B0_ra_resampled,-1);

	TH1D* h1_kick1_B2_B0_ra = (TH1D*)h1_kick1_B2_ra->Clone("h1_kick1_B2_B0_ra");
	h1_kick1_B2_B0_ra->Add(h1_kick1_B0_ra_resampled,-1);

	TH1D* h1_kick1_B3_B0_ra = (TH1D*)h1_kick1_B3_ra->Clone("h1_kick1_B3_B0_ra");
	h1_kick1_B3_B0_ra->Add(h1_kick1_B0_ra_resampled,-1);




	TH2D* h1_kick1_B1_B0_ra_spectrograph = getSpectrograph(h1_kick1_B1_B0_ra, 0, 8, 0.5);
	TH2D* h1_kick1_B1_0_ra_spectrograph = getSpectrograph(h1_kick1_B1_0_ra, 0, 8, 0.5);


//////////////////////// Style


	h1_kick1_B0_1_ra->SetTitle("Kick 1 (0 A) (Oct 14)");
	h1_kick1_B0_2_ra->SetTitle("Kick 1 (0 A) (Oct 23)");
	h1_kick1_B1_p_ra->SetTitle("Kick 1 (5173 A, Pos) (Jan 26)");
	h1_kick1_B1_n_ra->SetTitle("Kick 1 (5173 A, Neg) (Jan 28)");
	h1_kick1_B1_0_ra->SetTitle("Kick 1 (5173 A, Zero) (Jan 25)");

	h1_kick1_B0_ra->SetTitle("Kick 1 (0 A)");
	h1_kick1_B1_ra->SetTitle("Kick 1 (5173 A)");
	h1_kick1_B2_ra->SetTitle("Kick 1 (3043 A)");
	h1_kick1_B3_ra->SetTitle("Kick 1 (4353 A)");

	h1_kick1_B1_B0_ra->SetTitle("5173A - 0A");
	h1_kick1_B2_B0_ra->SetTitle("3043A - 0A");
	h1_kick1_B3_B0_ra->SetTitle("4353A - 0A");


	h1_kick1_B0_ra->SetLineWidth(2);
	h1_kick1_B1_ra->SetLineWidth(2);
	h1_kick1_B2_ra->SetLineWidth(2);
	h1_kick1_B3_ra->SetLineWidth(2);

	h1_kick1_B1_B0_ra->SetLineWidth(2);
	h1_kick1_B2_B0_ra->SetLineWidth(2);
	h1_kick1_B3_B0_ra->SetLineWidth(2);



	h1_kick1_B0_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B1_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B2_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B3_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B0_ra->GetYaxis()->SetRangeUser(-50,30);
	h1_kick1_B1_ra->GetYaxis()->SetRangeUser(-50,30);
	h1_kick1_B2_ra->GetYaxis()->SetRangeUser(-50,30);
	h1_kick1_B3_ra->GetYaxis()->SetRangeUser(-50,30);

	h1_kick1_B1_B0_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B2_B0_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B3_B0_ra->GetXaxis()->SetRangeUser(0,1.5);
	h1_kick1_B1_B0_ra->GetYaxis()->SetRangeUser(-50,30);
	h1_kick1_B2_B0_ra->GetYaxis()->SetRangeUser(-50,30);
	h1_kick1_B3_B0_ra->GetYaxis()->SetRangeUser(-50,30);


	h1_kick1_B0_1_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_2_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_p_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_n_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_1_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_2_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_p_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_n_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_1_ra->SetLineWidth(2);
	h1_kick1_B0_2_ra->SetLineWidth(2);
	h1_kick1_B1_p_ra->SetLineWidth(2);
	h1_kick1_B1_n_ra->SetLineWidth(2);
	h1_kick1_B0_1_ra->SetLineColor(kBlue);
	h1_kick1_B0_2_ra->SetLineColor(kRed);
	h1_kick1_B1_p_ra->SetLineColor(kBlue);
	h1_kick1_B1_n_ra->SetLineColor(kRed);



	h1_kick1_B1_0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_0_ra->SetLineWidth(2);
	h1_kick1_B1_0_ra->SetLineColor(kBlack);


	h1_kick1_B0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B0_ra->SetLineWidth(2);
	h1_kick1_B1_ra->SetLineWidth(2);
	h1_kick1_B0_ra->SetLineColor(kBlue);
	h1_kick1_B1_ra->SetLineColor(kRed);


	h1_kick1_B1_B0_ra->GetXaxis()->SetRangeUser(-1,1);
	h1_kick1_B1_B0_ra->GetYaxis()->SetRangeUser(-100,150);
	h1_kick1_B1_B0_ra->SetLineWidth(2);
	h1_kick1_B1_B0_ra->SetLineColor(kRed);




//////////////////////// Drawing


	gStyle->SetOptStat(0);

	new TCanvas();
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B0_2_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_1_ra->SetLineColor(1);
	h1_kick1_B0_2_ra->SetLineColor(2);
	h1_kick1_B0a_ra->SetLineColor(3);
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B0_2_ra->DrawCopy("HIST SAME");
	h1_kick1_B0a_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_B0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_B1_B0_ra->DrawCopy("HIST");
	h1_kick1_B1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	
	new TCanvas();
	h1_kick1_B0_ra_resampled->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas("","",1200,800);
	h1_kick1_B0_1_ra->SetLineWidth(2);
	h1_kick1_B0_1_ra->SetLineColor(kBlack);
	h1_kick1_B1_ra->SetLineWidth(2);
	h1_kick1_B1_ra->SetLineColor(kRed);
	h1_kick1_B0_1_ra->GetXaxis()->SetRangeUser(0,1);
	h1_kick1_B0_1_ra->GetYaxis()->SetRangeUser(-30,20);
	h1_kick1_B0_1_ra->DrawCopy("HIST");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
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
	h1_kick1_B1_0_ra_spectrograph->SetContour(100);
	h1_kick1_B1_0_ra_spectrograph->DrawCopy("colz");




	new TCanvas();
	h1_kick1_B0_ra->SetLineColor(1);
	h1_kick1_B2_ra->SetLineColor(2);
	h1_kick1_B3_ra->SetLineColor(3);
	h1_kick1_B1_ra->SetLineColor(4);
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B2_ra->DrawCopy("HIST SAME");
	h1_kick1_B3_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_B2_B0_ra->SetLineColor(2);
	h1_kick1_B3_B0_ra->SetLineColor(3);
	h1_kick1_B1_B0_ra->SetLineColor(4);
	h1_kick1_B2_B0_ra->DrawCopy("HIST");
	h1_kick1_B3_B0_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_B0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();







	///////////////// Fitting


	TF1* f_exp = new TF1("f_exp","[0]-[1]*exp(-x/[2])");
	f_exp->SetParameters(0,10,0.07);
	f_exp->SetParNames("offset","A","#tau");

	TF1* f_exp_17 = new TF1("f_exp_17","[0]-[1]*exp(-x/[2])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])");
	f_exp_17->SetParameters(0,10,0.07,2,0.3,17,0);
	f_exp_17->SetParLimits(3,0.5,10);
	f_exp_17->SetParLimits(4,0.1,1.0);
	f_exp_17->SetParLimits(5,15,19);
	f_exp_17->SetParNames("offset","A","#tau","A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	TF1* f_exp_17_osc = new TF1("f_exp_17_osc","[0]-[1]*exp(-x/[2])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])+[7]*exp(-x/[8])*sin([9]*6.283185307*x+[10])");
	f_exp_17_osc->SetParameters(0,10,0.07,2,0.3,17,0,5,0.5,3);
	f_exp_17_osc->SetParLimits(3,0.5,10);
	f_exp_17_osc->SetParLimits(4,0.1,1.0);
	f_exp_17_osc->SetParLimits(5,15,19);
	f_exp_17_osc->SetParLimits(7,0.1,20);
	f_exp_17_osc->SetParLimits(8,0.1,2);
	f_exp_17_osc->SetParLimits(9,1.5,5);
	f_exp_17_osc->SetParNames("offset","A","#tau","A_{1}","#tau_{1}","f_{1}","#phi_{1}","A_{2}","#tau_{2}","f_{2}");
	f_exp_17_osc->SetParName(10,"#phi_{2}");


	TGraph* g_chi2 = new TGraph();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	TGraphErrors* g_amp17 = new TGraphErrors();
	TGraphErrors* g_tau17 = new TGraphErrors();
	TGraphErrors* g_freq17 = new TGraphErrors();
	TGraphErrors* g_amp2 = new TGraphErrors();
	TGraphErrors* g_tau2 = new TGraphErrors();
	TGraphErrors* g_freq2 = new TGraphErrors();




	gStyle->SetOptFit(111);

	new TCanvas();
	h1_kick1_B0_ra->Draw("HIST");
	h1_kick1_B0_ra->Fit(f_exp,"S+","",0.02,0.7);
	f_exp->DrawCopy("SAME");
	gPad->SetGridy();


	for (int k=0; k<f_exp->GetNpar(); k++){
		f_exp_17->SetParameter(k,f_exp->GetParameter(k));
	}
	new TCanvas();
	h1_kick1_B0_ra->Draw("HIST");
	h1_kick1_B0_ra->Fit(f_exp_17,"S+","",0.02,0.7);
	f_exp_17->DrawCopy("SAME");
	gPad->SetGridy();

	for (int k=0; k<f_exp_17->GetNpar(); k++){
		f_exp_17_osc->SetParameter(k,f_exp_17->GetParameter(k));
	}
	new TCanvas();
	h1_kick1_B1_ra->Draw("HIST");
	TFitResultPtr h1_kick1_B1_ra_fit = h1_kick1_B1_ra->Fit(f_exp_17_osc,"S+","",0.02,0.7);
	f_exp_17_osc->DrawCopy("SAME");
	gPad->SetGridy();
	new TCanvas();
	h1_kick1_B2_ra->Draw("HIST");
	TFitResultPtr h1_kick1_B2_ra_fit = h1_kick1_B2_ra->Fit(f_exp_17_osc,"S+","",0.02,0.7);
	f_exp_17_osc->DrawCopy("SAME");
	gPad->SetGridy();
	new TCanvas();
	h1_kick1_B3_ra->Draw("HIST");
	TFitResultPtr h1_kick1_B3_ra_fit = h1_kick1_B3_ra->Fit(f_exp_17_osc,"S+","",0.02,0.7);
	f_exp_17_osc->DrawCopy("SAME");
	gPad->SetGridy();
	new TCanvas();
	h1_kick1_B0_ra->Draw("HIST");
	TFitResultPtr h1_kick1_B0_ra_fit = h1_kick1_B0_ra->Fit(f_exp_17,"S+","",0.02,0.7);
	f_exp_17->DrawCopy("SAME");
	gPad->SetGridy();


	g_chi2->SetPoint(0,3043,h1_kick1_B2_ra_fit->Chi2()/h1_kick1_B2_ra_fit->Ndf());
	g_amp->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(1));
	g_tau->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(2));
	g_amp17->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(3));
	g_tau17->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(4));
	g_freq17->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(5));
	g_amp2->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(7));
	g_tau2->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(8));
	g_freq2->SetPoint(0,3043,h1_kick1_B2_ra_fit->Parameter(9));
	g_amp->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(1));
	g_tau->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(2));
	g_amp17->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(3));
	g_tau17->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(4));
	g_freq17->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(5));
	g_amp2->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(7));
	g_tau2->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(8));
	g_freq2->SetPointError(0,0,h1_kick1_B2_ra_fit->ParError(9));

	g_chi2->SetPoint(1,4353,h1_kick1_B3_ra_fit->Chi2()/h1_kick1_B2_ra_fit->Ndf());
	g_amp->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(1));
	g_tau->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(2));
	g_amp17->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(3));
	g_tau17->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(4));
	g_freq17->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(5));
	g_amp2->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(7));
	g_tau2->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(8));
	g_freq2->SetPoint(1,4353,h1_kick1_B3_ra_fit->Parameter(9));
	g_amp->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(1));
	g_tau->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(2));
	g_amp17->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(3));
	g_tau17->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(4));
	g_freq17->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(5));
	g_amp2->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(7));
	g_tau2->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(8));
	g_freq2->SetPointError(1,0,h1_kick1_B3_ra_fit->ParError(9));

	g_chi2->SetPoint(2,5173,h1_kick1_B1_ra_fit->Chi2()/h1_kick1_B2_ra_fit->Ndf());
	g_amp->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(1));
	g_tau->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(2));
	g_amp17->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(3));
	g_tau17->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(4));
	g_freq17->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(5));
	g_amp2->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(7));
	g_tau2->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(8));
	g_freq2->SetPoint(2,5173,h1_kick1_B1_ra_fit->Parameter(9));
	g_amp->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(1));
	g_tau->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(2));
	g_amp17->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(3));
	g_tau17->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(4));
	g_freq17->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(5));
	g_amp2->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(7));
	g_tau2->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(8));
	g_freq2->SetPointError(2,0,h1_kick1_B1_ra_fit->ParError(9));

	g_chi2->SetPoint(3,0,h1_kick1_B0_ra_fit->Chi2()/h1_kick1_B2_ra_fit->Ndf());
	g_amp->SetPoint(3,0,h1_kick1_B0_ra_fit->Parameter(1));
	g_tau->SetPoint(3,0,h1_kick1_B0_ra_fit->Parameter(2));
	g_amp17->SetPoint(3,0,h1_kick1_B0_ra_fit->Parameter(3));
	g_tau17->SetPoint(3,0,h1_kick1_B0_ra_fit->Parameter(4));
	g_freq17->SetPoint(3,0,h1_kick1_B0_ra_fit->Parameter(5));
	g_amp->SetPointError(3,0,h1_kick1_B0_ra_fit->ParError(1));
	g_tau->SetPointError(3,0,h1_kick1_B0_ra_fit->ParError(2));
	g_amp17->SetPointError(3,0,h1_kick1_B0_ra_fit->ParError(3));
	g_tau17->SetPointError(3,0,h1_kick1_B0_ra_fit->ParError(4));
	g_freq17->SetPointError(3,0,h1_kick1_B0_ra_fit->ParError(5));


	new TCanvas();
	h1_kick1_B0_ra->SetLineColor(1);
	h1_kick1_B2_ra->SetLineColor(2);
	h1_kick1_B3_ra->SetLineColor(3);
	h1_kick1_B1_ra->SetLineColor(4);
	h1_kick1_B0_ra->DrawCopy("HIST");
	h1_kick1_B2_ra->DrawCopy("HIST SAME");
	h1_kick1_B3_ra->DrawCopy("HIST SAME");
	h1_kick1_B1_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();
	h1_kick1_B0_ra->GetFunction("f_exp_17")->DrawCopy("SAME");
	h1_kick1_B2_ra->GetFunction("f_exp_17_osc")->DrawCopy("SAME");
	h1_kick1_B3_ra->GetFunction("f_exp_17_osc")->DrawCopy("SAME");
	h1_kick1_B1_ra->GetFunction("f_exp_17_osc")->DrawCopy("SAME");



	g_chi2->Sort();
	g_amp->Sort();
	g_tau->Sort();
	g_amp17->Sort();
	g_tau17->Sort();
	g_freq17->Sort();
	g_amp2->Sort();
	g_tau2->Sort();
	g_freq2->Sort();
	g_chi2->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_amp17->SetMarkerStyle(20);
	g_tau17->SetMarkerStyle(20);
	g_freq17->SetMarkerStyle(20);
	g_amp2->SetMarkerStyle(20);
	g_tau2->SetMarkerStyle(20);
	g_freq2->SetMarkerStyle(20);
	g_chi2->GetXaxis()->SetTitle("Magnet current [A]");
	g_amp->GetXaxis()->SetTitle("Magnet current [A]");
	g_tau->GetXaxis()->SetTitle("Magnet current [A]");
	g_amp17->GetXaxis()->SetTitle("Magnet current [A]");
	g_tau17->GetXaxis()->SetTitle("Magnet current [A]");
	g_freq17->GetXaxis()->SetTitle("Magnet current [A]");
	g_amp2->GetXaxis()->SetTitle("Magnet current [A]");
	g_tau2->GetXaxis()->SetTitle("Magnet current [A]");
	g_freq2->GetXaxis()->SetTitle("Magnet current [A]");
	g_chi2->GetYaxis()->SetTitle("Chi2/ndf");
	g_amp->GetYaxis()->SetTitle("amp");
	g_tau->GetYaxis()->SetTitle("tau");
	g_amp17->GetYaxis()->SetTitle("amp17");
	g_tau17->GetYaxis()->SetTitle("tau17");
	g_freq17->GetYaxis()->SetTitle("freq17");
	g_amp2->GetYaxis()->SetTitle("amp2");
	g_tau2->GetYaxis()->SetTitle("tau2");
	g_freq2->GetYaxis()->SetTitle("freq2");


	new TCanvas();
	g_chi2->Draw("APLZ");
	new TCanvas();
	g_amp->Draw("APLZ");
	new TCanvas();
	g_tau->Draw("APLZ");
	new TCanvas();
	g_amp17->Draw("APLZ");
	new TCanvas();
	g_tau17->Draw("APLZ");
	new TCanvas();
	g_freq17->Draw("APLZ");
	new TCanvas();
	g_amp2->Draw("APLZ");
	new TCanvas();
	g_tau2->Draw("APLZ");
	new TCanvas();
	g_freq2->Draw("APLZ");

}
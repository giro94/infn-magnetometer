#include "../analysis_tools.C"

void plot_INFN_UMass_paper(){

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10;

	TFile* fin = TFile::Open("INFN_UMass.root");

	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");
	f_blumlein->SetParameters(125.0,-0.3,-1000.0);
	double xmin_blum = -0.42;
	double xmax_blum = -0.22;

	TH1D* h1_kick1_R0_ra = (TH1D*)fin->Get("h1_kick1_R0_ra");
	TH1D* h1_kick1_R1_ra = (TH1D*)fin->Get("h1_kick1_R1_ra");

	double blum_R0 = 126.397;
	double blum_R0err = 2.87448;
	double blum_R1 = 158.218;
	double blum_R1err = 3.59815;
	double rel_err_INFN_R0 = blum_R0err/blum_R0;
	double rel_err_INFN_R1 = blum_R1err/blum_R1;

	f_blumlein->SetParameters(125.0,-0.32,-3000.0);
	TFitResultPtr fit_blum_R0 = h1_kick1_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_R1 = h1_kick1_R1_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	double blum_fit_R0 = fit_blum_R0->Parameter(0);
	double blum_fit_R1 = fit_blum_R1->Parameter(0);
	h1_kick1_R0_ra->Scale(blum_R0/blum_fit_R0);
	h1_kick1_R1_ra->Scale(blum_R1/blum_fit_R1);

	TH1D* h1_K1_R0 = (TH1D*)fin->Get("h1_K1_R0");
	TH1D* h1_K1_R3p2 = (TH1D*)fin->Get("h1_K1_R3p2");
	TH1D* h1_K1_R6p6 = (TH1D*)fin->Get("h1_K1_R6p6");
	TH1D* h1_K1_Rm6p6 = (TH1D*)fin->Get("h1_K1_Rm6p6");
	TH1D* h1_K3_R0 = (TH1D*)fin->Get("h1_K3_R0");



	gROOT->SetBatch(kTRUE);


	TGraphErrors* g_INFN_blumlein = new TGraphErrors();
	TGraphErrors* g_UMass_blumlein = new TGraphErrors();
	TGraphErrors* g_all_blumlein = new TGraphErrors();
	TGraphErrors* g_INFN_amp30 = new TGraphErrors();
	TGraphErrors* g_UMass_amp30 = new TGraphErrors();
	TGraphErrors* g_all_amp30 = new TGraphErrors();
	TGraphErrors* g_INFN_ampexp = new TGraphErrors();
	TGraphErrors* g_UMass_ampexp = new TGraphErrors();
	TGraphErrors* g_all_ampexp = new TGraphErrors();
	TGraphErrors* g_INFN_ampintegral = new TGraphErrors();
	TGraphErrors* g_UMass_ampintegral = new TGraphErrors();
	TGraphErrors* g_all_ampintegral = new TGraphErrors();

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

	/////////Normalize according to 30/01/2025 email from David
	double blum_K1_R0 = 123.0;
	double blum_K1_R0err = 1.3;
	double blum_K1_R3p2 = 120.7;
	double blum_K1_R3p2err = 1.9;
	double blum_K1_R6p6 = 126.7;
	double blum_K1_R6p6err = 1.6;
	double blum_K1_Rm6p6 = 121.0;
	double blum_K1_Rm6p6err = 1.3;
	double blum_K3_R0 = 126.3;
	double blum_K3_R0err = 1.3;
	double rel_err_UMass_K1_R0 = blum_K1_R0err/blum_K1_R0;
	double rel_err_UMass_K1_R3p2 = blum_K1_R3p2err/blum_K1_R3p2;
	double rel_err_UMass_K1_R6p6 = blum_K1_R6p6err/blum_K1_R6p6;
	double rel_err_UMass_K1_Rm6p6 = blum_K1_Rm6p6err/blum_K1_Rm6p6;
	double rel_err_UMass_K3_R0 = blum_K3_R0err/blum_K3_R0;
	f_blumlein->SetParameters(125.0,-0.32,-3000.0);
	TFitResultPtr fit_blum_K1_R0 = h1_K1_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K1_R3p2 = h1_K1_R3p2_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K1_R6p6 = h1_K1_R6p6_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K1_Rm6p6 = h1_K1_Rm6p6_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_K3_R0 = h1_K3_R0_ra->Fit("f_blumlein","WS0N","",xmin_blum,xmax_blum);
	double blum_fit_K1_R0 = fit_blum_K1_R0->Parameter(0);
	double blum_fit_K1_R3p2 = fit_blum_K1_R3p2->Parameter(0);
	double blum_fit_K1_R6p6 = fit_blum_K1_R6p6->Parameter(0);
	double blum_fit_K1_Rm6p6 = fit_blum_K1_Rm6p6->Parameter(0);
	double blum_fit_K3_R0 = fit_blum_K3_R0->Parameter(0);
	h1_K1_R0_ra->Scale(blum_K1_R0/blum_fit_K1_R0);
	h1_K1_R3p2_ra->Scale(blum_K1_R3p2/blum_fit_K1_R3p2);
	h1_K1_R6p6_ra->Scale(blum_K1_R6p6/blum_fit_K1_R6p6);
	h1_K1_Rm6p6_ra->Scale(blum_K1_Rm6p6/blum_fit_K1_Rm6p6);
	h1_K3_R0_ra->Scale(blum_K3_R0/blum_fit_K3_R0);
	/////////


	//Rebin UMass hists to match INFN ones
	cout<<"\nRebin UMass hists\n";

	TH1D* h1_K1_R0_ra_rebinned = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R0_ra_rebinned");
	TH1D* h1_K1_R3p2_ra_rebinned = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R3p2_ra_rebinned");
	TH1D* h1_K1_R6p6_ra_rebinned = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_R6p6_ra_rebinned");
	TH1D* h1_K1_Rm6p6_ra_rebinned = (TH1D*)h1_kick1_R0_ra->Clone("h1_K1_Rm6p6_ra_rebinned");
	TH1D* h1_K3_R0_ra_rebinned = (TH1D*)h1_kick1_R0_ra->Clone("h1_K3_R0_ra_rebinned");

	h1_K1_R0_ra_rebinned->Reset();
	h1_K1_R3p2_ra_rebinned->Reset();
	h1_K1_R6p6_ra_rebinned->Reset();
	h1_K1_Rm6p6_ra_rebinned->Reset();
	h1_K3_R0_ra_rebinned->Reset();
	for (int bn=1; bn<=h1_kick1_R0_ra->GetNbinsX(); bn++){
		double x = h1_kick1_R0_ra->GetBinCenter(bn);
		h1_K1_R0_ra_rebinned->Fill(x,h1_K1_R0_ra->Interpolate(x));
		h1_K1_R3p2_ra_rebinned->Fill(x,h1_K1_R3p2_ra->Interpolate(x));
		h1_K1_R6p6_ra_rebinned->Fill(x,h1_K1_R6p6_ra->Interpolate(x));
		h1_K1_Rm6p6_ra_rebinned->Fill(x,h1_K1_Rm6p6_ra->Interpolate(x));
		h1_K3_R0_ra_rebinned->Fill(x,h1_K3_R0_ra->Interpolate(x));
	}
	h1_K1_R0_ra_rebinned->SetTitle("UMass K1 0.0 mm");
	h1_K1_R3p2_ra_rebinned->SetTitle("UMass K1 3.2 mm");
	h1_K1_R6p6_ra_rebinned->SetTitle("UMass K1 6.6 mm");
	h1_K1_Rm6p6_ra_rebinned->SetTitle("UMass K1 -6.6 mm");
	h1_K3_R0_ra_rebinned->SetTitle("UMass K3 0.0 mm");

	TH1D* h1_K1_R0_ra_rescaled = (TH1D*)h1_K1_R0_ra_rebinned->Clone("h1_K1_R0_ra_rescaled");
	TH1D* h1_K1_R3p2_ra_rescaled = (TH1D*)h1_K1_R3p2_ra_rebinned->Clone("h1_K1_R3p2_ra_rescaled");
	TH1D* h1_K1_R6p6_ra_rescaled = (TH1D*)h1_K1_R6p6_ra_rebinned->Clone("h1_K1_R6p6_ra_rescaled");
	TH1D* h1_K1_Rm6p6_ra_rescaled = (TH1D*)h1_K1_Rm6p6_ra_rebinned->Clone("h1_K1_Rm6p6_ra_rescaled");


	//Rescale K1 to match K3 strength
	cout<<"\nRescake K1\n";

	double k1_to_k3 = 55.0/53.1;

	h1_K1_R0_ra_rescaled->Scale(k1_to_k3);
	h1_K1_R3p2_ra_rescaled->Scale(k1_to_k3);
	h1_K1_R6p6_ra_rescaled->Scale(k1_to_k3);
	h1_K1_Rm6p6_ra_rescaled->Scale(k1_to_k3);

	double dx_INFN = 2;
	double dx_UMass = 3;

	//Fill the graphs

	//blumlein
	cout<<"\nFitting blumleins\n";

	f_blumlein->SetParameters(125.0,-0.32,-3000.0);
	TFitResultPtr fit_blum_INFN_R0 = h1_kick1_R0_ra->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_INFN_R1 = h1_kick1_R1_ra->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_UMass_K1_R0 = h1_K1_R0_ra_rescaled->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_UMass_K1_R3p2 = h1_K1_R3p2_ra_rescaled->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_UMass_K1_R6p6 = h1_K1_R6p6_ra_rescaled->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_UMass_K1_Rm6p6 = h1_K1_Rm6p6_ra_rescaled->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);
	TFitResultPtr fit_blum_UMass_K3_R0 = h1_K3_R0_ra_rebinned->Fit("f_blumlein","WS","",xmin_blum,xmax_blum);

	double blum_INFN_R0 = fit_blum_INFN_R0->Parameter(0);
	double blum_INFN_R1 = fit_blum_INFN_R1->Parameter(0);
	double blum_UMass_K1_R0 = fit_blum_UMass_K1_R0->Parameter(0);
	double blum_UMass_K1_R3p2 = fit_blum_UMass_K1_R3p2->Parameter(0);
	double blum_UMass_K1_R6p6 = fit_blum_UMass_K1_R6p6->Parameter(0);
	double blum_UMass_K1_Rm6p6 = fit_blum_UMass_K1_Rm6p6->Parameter(0);
	double blum_UMass_K3_R0 = fit_blum_UMass_K3_R0->Parameter(0);

	double blumerr_INFN_R0 = fit_blum_INFN_R0->ParError(0);
	double blumerr_INFN_R1 = fit_blum_INFN_R1->ParError(0);
	double blumerr_UMass_K1_R0 = fit_blum_UMass_K1_R0->ParError(0);
	double blumerr_UMass_K1_R3p2 = fit_blum_UMass_K1_R3p2->ParError(0);
	double blumerr_UMass_K1_R6p6 = fit_blum_UMass_K1_R6p6->ParError(0);
	double blumerr_UMass_K1_Rm6p6 = fit_blum_UMass_K1_Rm6p6->ParError(0);
	double blumerr_UMass_K3_R0 = fit_blum_UMass_K3_R0->ParError(0);

	double calibErr_blum_INFN_R0 = blum_INFN_R0*rel_err_INFN_R0;
	double calibErr_blum_INFN_R1 = blum_INFN_R1*rel_err_INFN_R1;
	double calibErr_blum_UMassK1_0 = blum_UMass_K1_R0*rel_err_UMass_K1_R0;
	double calibErr_blum_UMassK1_3p2 = blum_UMass_K1_R3p2*rel_err_UMass_K1_R3p2;
	double calibErr_blum_UMassK1_6p6 = blum_UMass_K1_R6p6*rel_err_UMass_K1_R6p6;
	double calibErr_blum_UMassK1_m6p6 = blum_UMass_K1_Rm6p6*rel_err_UMass_K1_Rm6p6;
	double calibErr_blum_UMassK3_0 = blum_UMass_K3_R0*rel_err_UMass_K3_R0;
	
	g_INFN_blumlein->AddPoint(0,blum_INFN_R0);
	g_INFN_blumlein->SetPointError(g_INFN_blumlein->GetN()-1,dx_INFN,sqrt(calibErr_blum_INFN_R0*calibErr_blum_INFN_R0 + blumerr_INFN_R0*blumerr_INFN_R0));
	g_INFN_blumlein->AddPoint(17.5,blum_INFN_R1);
	g_INFN_blumlein->SetPointError(g_INFN_blumlein->GetN()-1,dx_INFN,sqrt(calibErr_blum_INFN_R1*calibErr_blum_INFN_R1 + blumerr_INFN_R1*blumerr_INFN_R1));
	
	g_UMass_blumlein->AddPoint(0,blum_UMass_K1_R0);
	g_UMass_blumlein->SetPointError(g_UMass_blumlein->GetN()-1,dx_UMass,sqrt(calibErr_blum_UMassK1_0*calibErr_blum_UMassK1_0 + blumerr_UMass_K1_R0*blumerr_UMass_K1_R0));
	g_UMass_blumlein->AddPoint(3.2,blum_UMass_K1_R3p2);
	g_UMass_blumlein->SetPointError(g_UMass_blumlein->GetN()-1,dx_UMass,sqrt(calibErr_blum_UMassK1_3p2*calibErr_blum_UMassK1_3p2 + blumerr_UMass_K1_R3p2*blumerr_UMass_K1_R3p2));
	g_UMass_blumlein->AddPoint(6.6,blum_UMass_K1_R6p6);
	g_UMass_blumlein->SetPointError(g_UMass_blumlein->GetN()-1,dx_UMass,sqrt(calibErr_blum_UMassK1_6p6*calibErr_blum_UMassK1_6p6 + blumerr_UMass_K1_R6p6*blumerr_UMass_K1_R6p6));
	g_UMass_blumlein->AddPoint(-6.6,blum_UMass_K1_Rm6p6);
	g_UMass_blumlein->SetPointError(g_UMass_blumlein->GetN()-1,dx_UMass,sqrt(calibErr_blum_UMassK1_m6p6*calibErr_blum_UMassK1_m6p6 + blumerr_UMass_K1_Rm6p6*blumerr_UMass_K1_Rm6p6));
	g_UMass_blumlein->AddPoint(0,blum_UMass_K3_R0);
	g_UMass_blumlein->SetPointError(g_UMass_blumlein->GetN()-1,dx_UMass,sqrt(calibErr_blum_UMassK3_0*calibErr_blum_UMassK3_0 + blumerr_UMass_K3_R0*blumerr_UMass_K3_R0));

	TList blum_list;
	blum_list.Add(g_INFN_blumlein);
	blum_list.Add(g_UMass_blumlein);
	g_all_blumlein->Merge(&blum_list);


	//transient @30µs
	cout<<"\nTransient at 30 us \n";

	double value_30_INFN_R0 = h1_kick1_R0_ra->Interpolate(0.03);
	double value_30_INFN_R1 = h1_kick1_R1_ra->Interpolate(0.03);
	double value_30_UMassK1_0 = h1_K1_R0_ra_rescaled->Interpolate(0.03);
	double value_30_UMassK1_3p2 = h1_K1_R3p2_ra_rescaled->Interpolate(0.03);
	double value_30_UMassK1_6p6 = h1_K1_R6p6_ra_rescaled->Interpolate(0.03);
	double value_30_UMassK1_m6p6 = h1_K1_Rm6p6_ra_rescaled->Interpolate(0.03);
	double value_30_UMassK3_0 = h1_K3_R0_ra_rebinned->Interpolate(0.03);

	g_INFN_amp30->AddPoint(0,value_30_INFN_R0);
	g_INFN_amp30->SetPointError(g_INFN_amp30->GetN()-1,dx_INFN,abs(value_30_INFN_R0*rel_err_INFN_R0));
	g_INFN_amp30->AddPoint(17.5,value_30_INFN_R1);
	g_INFN_amp30->SetPointError(g_INFN_amp30->GetN()-1,dx_INFN,abs(value_30_INFN_R1*rel_err_INFN_R1));

	g_UMass_amp30->AddPoint(0,value_30_UMassK1_0);
	g_UMass_amp30->SetPointError(g_UMass_amp30->GetN()-1,dx_UMass,abs(value_30_UMassK1_0*rel_err_UMass_K1_R0));
	g_UMass_amp30->AddPoint(3.2,value_30_UMassK1_3p2);
	g_UMass_amp30->SetPointError(g_UMass_amp30->GetN()-1,dx_UMass,abs(value_30_UMassK1_3p2*rel_err_UMass_K1_R3p2));
	g_UMass_amp30->AddPoint(6.6,value_30_UMassK1_6p6);
	g_UMass_amp30->SetPointError(g_UMass_amp30->GetN()-1,dx_UMass,abs(value_30_UMassK1_6p6*rel_err_UMass_K1_R6p6));
	g_UMass_amp30->AddPoint(-6.6,value_30_UMassK1_m6p6);
	g_UMass_amp30->SetPointError(g_UMass_amp30->GetN()-1,dx_UMass,abs(value_30_UMassK1_m6p6*rel_err_UMass_K1_Rm6p6));
	g_UMass_amp30->AddPoint(0,value_30_UMassK3_0);
	g_UMass_amp30->SetPointError(g_UMass_amp30->GetN()-1,dx_UMass,abs(value_30_UMassK3_0*rel_err_UMass_K3_R0));

	TList amp30_list;
	amp30_list.Add(g_INFN_amp30);
	amp30_list.Add(g_UMass_amp30);
	g_all_amp30->Merge(&amp30_list);

	//transient exp amp
	cout<<"\nFitting transient amplitude \n";

	double xmin_exp = 0.015;
	double xmax_exp = 0.7;

	TF1* f_exp = new TF1("f_exp","[2]+[0]*exp(-(x-0.03)/[1])",xmin_exp,xmax_exp);
	f_exp->SetParameters(-15,0.05,0);
	f_exp->FixParameter(2,0);
	f_exp->SetParNames("A30","#tau","offset");
	TF1* f_exp_17 = new TF1("f_exp_17","[2]+[0]*exp(-(x-0.03)/[1])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(-15,0.05,0,2,0.3,17,0);
	f_exp_17->FixParameter(2,0);
	f_exp_17->SetParLimits(3,0.5,10);
	f_exp_17->SetParLimits(4,0.1,1.0);
	f_exp_17->SetParLimits(5,15,19);
	f_exp_17->SetParNames("A30","#tau","offset","A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	new TCanvas();
	h1_kick1_R0_ra->GetXaxis()->SetRangeUser(0,1);
	h1_kick1_R0_ra->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_INFN_R0 = h1_kick1_R0_ra->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_INFN_R0 = h1_kick1_R0_ra->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_kick1_R1_ra->GetXaxis()->SetRangeUser(0,1);
	h1_kick1_R1_ra->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_INFN_R1 = h1_kick1_R1_ra->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_INFN_R1 = h1_kick1_R1_ra->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_K1_R0_ra_rescaled->GetXaxis()->SetRangeUser(0,1);
	h1_K1_R0_ra_rescaled->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_UMassK1_0 = h1_K1_R0_ra_rescaled->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_UMassK1_0 = h1_K1_R0_ra_rescaled->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_K1_R3p2_ra_rescaled->GetXaxis()->SetRangeUser(0,1);
	h1_K1_R3p2_ra_rescaled->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_UMassK1_3p2 = h1_K1_R3p2_ra_rescaled->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_UMassK1_3p2 = h1_K1_R3p2_ra_rescaled->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_K1_R6p6_ra_rescaled->GetXaxis()->SetRangeUser(0,1);
	h1_K1_R6p6_ra_rescaled->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_UMassK1_6p6 = h1_K1_R6p6_ra_rescaled->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_UMassK1_6p6 = h1_K1_R6p6_ra_rescaled->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_K1_Rm6p6_ra_rescaled->GetXaxis()->SetRangeUser(0,1);
	h1_K1_Rm6p6_ra_rescaled->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_UMassK1_m6p6 = h1_K1_Rm6p6_ra_rescaled->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_UMassK1_m6p6 = h1_K1_Rm6p6_ra_rescaled->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");
	new TCanvas();
	h1_K3_R0_ra_rebinned->GetXaxis()->SetRangeUser(0,1);
	h1_K3_R0_ra_rebinned->Draw("HIST");
	f_exp->SetParameters(-15,0.05,0);
	TFitResultPtr fit_exp_UMassK3_0 = h1_K3_R0_ra_rebinned->Fit("f_exp","WS","",xmin_exp,xmax_exp);
	//f_exp_17->SetParameters(15,0.05,0,2,0.3,17,0);
	//TFitResultPtr fit_exp_17_UMassK3_0 = h1_K3_R0_ra_rebinned->Fit("f_exp_17","WS","",xmin_exp,xmax_exp);
	f_exp->DrawCopy("SAME");
	//f_exp_17->DrawCopy("SAME");

	double value_exp_INFN_R0 = fit_exp_INFN_R0->Parameter(0);
	double value_exp_INFN_R1 = fit_exp_INFN_R1->Parameter(0);
	double value_exp_UMassK1_0 = fit_exp_UMassK1_0->Parameter(0);
	double value_exp_UMassK1_3p2 = fit_exp_UMassK1_3p2->Parameter(0);
	double value_exp_UMassK1_6p6 = fit_exp_UMassK1_6p6->Parameter(0);
	double value_exp_UMassK1_m6p6 = fit_exp_UMassK1_m6p6->Parameter(0);
	double value_exp_UMassK3_0 = fit_exp_UMassK3_0->Parameter(0);

	double valueErr_exp_INFN_R0 = fit_exp_INFN_R0->ParError(0);
	double valueErr_exp_INFN_R1 = fit_exp_INFN_R1->ParError(0);
	double valueErr_exp_UMassK1_0 = fit_exp_UMassK1_0->ParError(0);
	double valueErr_exp_UMassK1_3p2 = fit_exp_UMassK1_3p2->ParError(0);
	double valueErr_exp_UMassK1_6p6 = fit_exp_UMassK1_6p6->ParError(0);
	double valueErr_exp_UMassK1_m6p6 = fit_exp_UMassK1_m6p6->ParError(0);
	double valueErr_exp_UMassK3_0 = fit_exp_UMassK3_0->ParError(0);

	double calibErr_exp_INFN_R0 = value_exp_INFN_R0*rel_err_INFN_R0;
	double calibErr_exp_INFN_R1 = value_exp_INFN_R1*rel_err_INFN_R1;
	double calibErr_exp_UMassK1_0 = value_exp_UMassK1_0*rel_err_UMass_K1_R0;
	double calibErr_exp_UMassK1_3p2 = value_exp_UMassK1_3p2*rel_err_UMass_K1_R3p2;
	double calibErr_exp_UMassK1_6p6 = value_exp_UMassK1_6p6*rel_err_UMass_K1_R6p6;
	double calibErr_exp_UMassK1_m6p6 = value_exp_UMassK1_m6p6*rel_err_UMass_K1_Rm6p6;
	double calibErr_exp_UMassK3_0 = value_exp_UMassK3_0*rel_err_UMass_K3_R0;

	g_INFN_ampexp->AddPoint(0,value_exp_INFN_R0);
	g_INFN_ampexp->SetPointError(g_INFN_ampexp->GetN()-1,dx_INFN,sqrt(valueErr_exp_INFN_R0*valueErr_exp_INFN_R0 + calibErr_exp_INFN_R0*calibErr_exp_INFN_R0));
	g_INFN_ampexp->AddPoint(17.5,value_exp_INFN_R1);
	g_INFN_ampexp->SetPointError(g_INFN_ampexp->GetN()-1,dx_INFN,sqrt(valueErr_exp_INFN_R1*valueErr_exp_INFN_R1 + calibErr_exp_INFN_R1*calibErr_exp_INFN_R1));

	g_UMass_ampexp->AddPoint(0,value_exp_UMassK1_0);
	g_UMass_ampexp->SetPointError(g_UMass_ampexp->GetN()-1,dx_UMass,sqrt(valueErr_exp_UMassK1_0*valueErr_exp_UMassK1_0 + calibErr_exp_UMassK1_0*calibErr_exp_UMassK1_0));
	g_UMass_ampexp->AddPoint(3.2,value_exp_UMassK1_3p2);
	g_UMass_ampexp->SetPointError(g_UMass_ampexp->GetN()-1,dx_UMass,sqrt(valueErr_exp_UMassK1_3p2*valueErr_exp_UMassK1_3p2 + calibErr_exp_UMassK1_3p2*calibErr_exp_UMassK1_3p2));
	g_UMass_ampexp->AddPoint(6.6,value_exp_UMassK1_6p6);
	g_UMass_ampexp->SetPointError(g_UMass_ampexp->GetN()-1,dx_UMass,sqrt(valueErr_exp_UMassK1_6p6*valueErr_exp_UMassK1_6p6 + calibErr_exp_UMassK1_6p6*calibErr_exp_UMassK1_6p6));
	g_UMass_ampexp->AddPoint(-6.6,value_exp_UMassK1_m6p6);
	g_UMass_ampexp->SetPointError(g_UMass_ampexp->GetN()-1,dx_UMass,sqrt(valueErr_exp_UMassK1_m6p6*valueErr_exp_UMassK1_m6p6 + calibErr_exp_UMassK1_m6p6*calibErr_exp_UMassK1_m6p6));
	g_UMass_ampexp->AddPoint(0,value_exp_UMassK3_0);
	g_UMass_ampexp->SetPointError(g_UMass_ampexp->GetN()-1,dx_UMass,sqrt(valueErr_exp_UMassK3_0*valueErr_exp_UMassK3_0 + calibErr_exp_UMassK3_0*calibErr_exp_UMassK3_0));

	TList ampexp_list;
	ampexp_list.Add(g_INFN_ampexp);
	ampexp_list.Add(g_UMass_ampexp);
	g_all_ampexp->Merge(&ampexp_list);

	cout<<"value_exp_INFN_R0: " << value_exp_INFN_R0 << " mG\n";
	cout<<"value_exp_INFN_R1: " << value_exp_INFN_R1 << " mG\n";
	cout<<"value_exp_UMassK1_0: " << value_exp_UMassK1_0 << " mG\n";
	cout<<"value_exp_UMassK1_3p2: " << value_exp_UMassK1_3p2 << " mG\n";
	cout<<"value_exp_UMassK1_6p6: " << value_exp_UMassK1_6p6 << " mG\n";
	cout<<"value_exp_UMassK1_m6p6: " << value_exp_UMassK1_m6p6 << " mG\n";
	cout<<"value_exp_UMassK3_0: " << value_exp_UMassK3_0 << " mG\n";



	//transient integral
	cout<<"\nTransient integral \n";

	int firstbin = h1_kick1_R0_ra->FindBin(0.03);
	int lastbin = h1_kick1_R0_ra->FindBin(0.7);

	double value_integral_INFN_R0 = h1_kick1_R0_ra->Integral(firstbin,lastbin);
	double value_integral_INFN_R1 = h1_kick1_R1_ra->Integral(firstbin,lastbin);
	double value_integral_UMassK1_0 = h1_K1_R0_ra_rescaled->Integral(firstbin,lastbin);
	double value_integral_UMassK1_3p2 = h1_K1_R3p2_ra_rescaled->Integral(firstbin,lastbin);
	double value_integral_UMassK1_6p6 = h1_K1_R6p6_ra_rescaled->Integral(firstbin,lastbin);
	double value_integral_UMassK1_m6p6 = h1_K1_Rm6p6_ra_rescaled->Integral(firstbin,lastbin);
	double value_integral_UMassK3_0 = h1_K3_R0_ra_rebinned->Integral(firstbin,lastbin);

	g_INFN_ampintegral->AddPoint(0,value_integral_INFN_R0);
	g_INFN_ampintegral->SetPointError(g_INFN_ampintegral->GetN()-1,dx_INFN,value_integral_INFN_R0*rel_err_INFN_R0);
	g_INFN_ampintegral->AddPoint(17.5,value_integral_INFN_R1);
	g_INFN_ampintegral->SetPointError(g_INFN_ampintegral->GetN()-1,dx_INFN,value_integral_INFN_R1*rel_err_INFN_R1);

	g_UMass_ampintegral->AddPoint(0,value_integral_UMassK1_0);
	g_UMass_ampintegral->SetPointError(g_UMass_ampintegral->GetN()-1,dx_UMass,value_integral_UMassK1_0*rel_err_UMass_K1_R0);
	g_UMass_ampintegral->AddPoint(3.2,value_integral_UMassK1_3p2);
	g_UMass_ampintegral->SetPointError(g_UMass_ampintegral->GetN()-1,dx_UMass,value_integral_UMassK1_3p2*rel_err_UMass_K1_R3p2);
	g_UMass_ampintegral->AddPoint(6.6,value_integral_UMassK1_6p6);
	g_UMass_ampintegral->SetPointError(g_UMass_ampintegral->GetN()-1,dx_UMass,value_integral_UMassK1_6p6*rel_err_UMass_K1_R6p6);
	g_UMass_ampintegral->AddPoint(-6.6,value_integral_UMassK1_m6p6);
	g_UMass_ampintegral->SetPointError(g_UMass_ampintegral->GetN()-1,dx_UMass,value_integral_UMassK1_m6p6*rel_err_UMass_K1_Rm6p6);
	g_UMass_ampintegral->AddPoint(0,value_integral_UMassK3_0);
	g_UMass_ampintegral->SetPointError(g_UMass_ampintegral->GetN()-1,dx_UMass,value_integral_UMassK3_0*rel_err_UMass_K3_R0);

	TList ampintegral_list;
	ampintegral_list.Add(g_INFN_ampintegral);
	ampintegral_list.Add(g_UMass_ampintegral);
	g_all_ampintegral->Merge(&ampintegral_list);



	//Normalized plots
	cout<<"Normalizing graphs\n";

	TGraphErrors* g_INFN_blumlein_normalized = (TGraphErrors*)g_INFN_blumlein->Clone("g_INFN_blumlein_normalized");
	TGraphErrors* g_UMass_blumlein_normalized = (TGraphErrors*)g_UMass_blumlein->Clone("g_UMass_blumlein_normalized");
	TGraphErrors* g_all_blumlein_normalized = new TGraphErrors();
	TGraphErrors* g_INFN_amp30_normalized = (TGraphErrors*)g_INFN_amp30->Clone("g_INFN_amp30_normalized");
	TGraphErrors* g_UMass_amp30_normalized = (TGraphErrors*)g_UMass_amp30->Clone("g_UMass_amp30_normalized");
	TGraphErrors* g_all_amp30_normalized = new TGraphErrors();
	TGraphErrors* g_INFN_ampexp_normalized = (TGraphErrors*)g_INFN_ampexp->Clone("g_INFN_ampexp_normalized");
	TGraphErrors* g_UMass_ampexp_normalized = (TGraphErrors*)g_UMass_ampexp->Clone("g_UMass_ampexp_normalized");
	TGraphErrors* g_all_ampexp_normalized = new TGraphErrors();
	TGraphErrors* g_INFN_ampintegral_normalized = (TGraphErrors*)g_INFN_ampintegral->Clone("g_INFN_ampintegral_normalized");
	TGraphErrors* g_UMass_ampintegral_normalized = (TGraphErrors*)g_UMass_ampintegral->Clone("g_UMass_ampintegral_normalized");
	TGraphErrors* g_all_ampintegral_normalized = new TGraphErrors();

	g_INFN_blumlein_normalized->Scale(1./g_INFN_blumlein->GetPointY(0));
	g_UMass_blumlein_normalized->Scale(1./g_UMass_blumlein->GetPointY(0));
	g_UMass_blumlein_normalized->RemovePoint(g_UMass_blumlein_normalized->GetN()-1); //Remove K3 point
	//Stupid root for scaling errors into negative values
	for (int i=0; i<g_INFN_blumlein_normalized->GetN(); i++) {g_INFN_blumlein_normalized->SetPointError(i,g_INFN_blumlein_normalized->GetErrorX(i),abs(g_INFN_blumlein_normalized->GetErrorY(i)));}
	for (int i=0; i<g_UMass_blumlein_normalized->GetN(); i++) {g_UMass_blumlein_normalized->SetPointError(i,g_UMass_blumlein_normalized->GetErrorX(i),abs(g_UMass_blumlein_normalized->GetErrorY(i)));}
	TList blum_norm_list;
	blum_norm_list.Add(g_INFN_blumlein_normalized);
	blum_norm_list.Add(g_UMass_blumlein_normalized);
	g_all_blumlein_normalized->Merge(&blum_norm_list);

	g_INFN_amp30_normalized->Scale(1./g_INFN_amp30->GetPointY(0));
	g_UMass_amp30_normalized->Scale(1./g_UMass_amp30->GetPointY(0));
	g_UMass_amp30_normalized->RemovePoint(g_UMass_amp30_normalized->GetN()-1); //Remove K3 point
	for (int i=0; i<g_INFN_amp30_normalized->GetN(); i++) {g_INFN_amp30_normalized->SetPointError(i,g_INFN_amp30_normalized->GetErrorX(i),abs(g_INFN_amp30_normalized->GetErrorY(i)));}
	for (int i=0; i<g_UMass_amp30_normalized->GetN(); i++) {g_UMass_amp30_normalized->SetPointError(i,g_UMass_amp30_normalized->GetErrorX(i),abs(g_UMass_amp30_normalized->GetErrorY(i)));}
	TList amp30_norm_list;
	amp30_norm_list.Add(g_INFN_amp30_normalized);
	amp30_norm_list.Add(g_UMass_amp30_normalized);
	g_all_amp30_normalized->Merge(&amp30_norm_list);

	g_INFN_ampexp_normalized->Scale(1./g_INFN_ampexp->GetPointY(0));
	g_UMass_ampexp_normalized->Scale(1./g_UMass_ampexp->GetPointY(0));
	g_UMass_ampexp_normalized->RemovePoint(g_UMass_ampexp_normalized->GetN()-1); //Remove K3 point
	for (int i=0; i<g_INFN_ampexp_normalized->GetN(); i++) {g_INFN_ampexp_normalized->SetPointError(i,g_INFN_ampexp_normalized->GetErrorX(i),abs(g_INFN_ampexp_normalized->GetErrorY(i)));}
	for (int i=0; i<g_UMass_ampexp_normalized->GetN(); i++) {g_UMass_ampexp_normalized->SetPointError(i,g_UMass_ampexp_normalized->GetErrorX(i),abs(g_UMass_ampexp_normalized->GetErrorY(i)));}
	TList ampexp_norm_list;
	ampexp_norm_list.Add(g_INFN_ampexp_normalized);
	ampexp_norm_list.Add(g_UMass_ampexp_normalized);
	g_all_ampexp_normalized->Merge(&ampexp_norm_list);

	g_INFN_ampintegral_normalized->Scale(1./g_INFN_ampintegral->GetPointY(0));
	g_UMass_ampintegral_normalized->Scale(1./g_UMass_ampintegral->GetPointY(0));
	g_UMass_ampintegral_normalized->RemovePoint(g_UMass_ampintegral_normalized->GetN()-1); //Remove K3 point
	for (int i=0; i<g_INFN_ampintegral_normalized->GetN(); i++) {g_INFN_ampintegral_normalized->SetPointError(i,g_INFN_ampintegral_normalized->GetErrorX(i),abs(g_INFN_ampintegral_normalized->GetErrorY(i)));}
	for (int i=0; i<g_UMass_ampintegral_normalized->GetN(); i++) {g_UMass_ampintegral_normalized->SetPointError(i,g_UMass_ampintegral_normalized->GetErrorX(i),abs(g_UMass_ampintegral_normalized->GetErrorY(i)));}
	TList ampintegral_norm_list;
	ampintegral_norm_list.Add(g_INFN_ampintegral_normalized);
	ampintegral_norm_list.Add(g_UMass_ampintegral_normalized);
	g_all_ampintegral_normalized->Merge(&ampintegral_norm_list);



	///////////////////////////////////


	g_INFN_blumlein->SetMarkerStyle(20);
	g_UMass_blumlein->SetMarkerStyle(20);
	g_all_blumlein->SetMarkerStyle(20);
	g_INFN_amp30->SetMarkerStyle(20);
	g_UMass_amp30->SetMarkerStyle(20);
	g_all_amp30->SetMarkerStyle(20);
	g_INFN_ampexp->SetMarkerStyle(20);
	g_UMass_ampexp->SetMarkerStyle(20);
	g_all_ampexp->SetMarkerStyle(20);
	g_INFN_ampintegral->SetMarkerStyle(20);
	g_UMass_ampintegral->SetMarkerStyle(20);
	g_all_ampintegral->SetMarkerStyle(20);

	g_INFN_blumlein_normalized->SetMarkerStyle(20);
	g_UMass_blumlein_normalized->SetMarkerStyle(20);
	g_all_blumlein_normalized->SetMarkerStyle(20);
	g_INFN_amp30_normalized->SetMarkerStyle(20);
	g_UMass_amp30_normalized->SetMarkerStyle(20);
	g_all_amp30_normalized->SetMarkerStyle(20);
	g_INFN_ampexp_normalized->SetMarkerStyle(20);
	g_UMass_ampexp_normalized->SetMarkerStyle(20);
	g_all_ampexp_normalized->SetMarkerStyle(20);
	g_INFN_ampintegral_normalized->SetMarkerStyle(20);
	g_UMass_ampintegral_normalized->SetMarkerStyle(20);
	g_all_ampintegral_normalized->SetMarkerStyle(20);

	g_INFN_blumlein->SetMarkerColor(kBlue);
	g_UMass_blumlein->SetMarkerColor(kRed);
	g_all_blumlein->SetMarkerColor(kBlack);
	g_INFN_amp30->SetMarkerColor(kBlue);
	g_UMass_amp30->SetMarkerColor(kRed);
	g_all_amp30->SetMarkerColor(kBlack);
	g_INFN_ampexp->SetMarkerColor(kBlue);
	g_UMass_ampexp->SetMarkerColor(kRed);
	g_all_ampexp->SetMarkerColor(kBlack);
	g_INFN_ampintegral->SetMarkerColor(kBlue);
	g_UMass_ampintegral->SetMarkerColor(kRed);
	g_all_ampintegral->SetMarkerColor(kBlack);

	g_INFN_blumlein_normalized->SetMarkerColor(kBlue);
	g_UMass_blumlein_normalized->SetMarkerColor(kRed);
	g_all_blumlein_normalized->SetMarkerColor(kBlack);
	g_INFN_amp30_normalized->SetMarkerColor(kBlue);
	g_UMass_amp30_normalized->SetMarkerColor(kRed);
	g_all_amp30_normalized->SetMarkerColor(kBlack);
	g_INFN_ampexp_normalized->SetMarkerColor(kBlue);
	g_UMass_ampexp_normalized->SetMarkerColor(kRed);
	g_all_ampexp_normalized->SetMarkerColor(kBlack);
	g_INFN_ampintegral_normalized->SetMarkerColor(kBlue);
	g_UMass_ampintegral_normalized->SetMarkerColor(kRed);
	g_all_ampintegral_normalized->SetMarkerColor(kBlack);

	g_INFN_blumlein->SetTitle("INFN");
	g_UMass_blumlein->SetTitle("UMass");
	g_all_blumlein->SetTitle("");
	g_INFN_amp30->SetTitle("INFN");
	g_UMass_amp30->SetTitle("UMass");
	g_all_amp30->SetTitle("");
	g_INFN_ampexp->SetTitle("INFN");
	g_UMass_ampexp->SetTitle("UMass");
	g_all_ampexp->SetTitle("");
	g_INFN_ampintegral->SetTitle("INFN");
	g_UMass_ampintegral->SetTitle("UMass");
	g_all_ampintegral->SetTitle("");

	g_INFN_blumlein_normalized->SetTitle("INFN");
	g_UMass_blumlein_normalized->SetTitle("UMass");
	g_all_blumlein_normalized->SetTitle("");
	g_INFN_amp30_normalized->SetTitle("INFN");
	g_UMass_amp30_normalized->SetTitle("UMass");
	g_all_amp30_normalized->SetTitle("");
	g_INFN_ampexp_normalized->SetTitle("INFN");
	g_UMass_ampexp_normalized->SetTitle("UMass");
	g_all_ampexp_normalized->SetTitle("");
	g_INFN_ampintegral_normalized->SetTitle("INFN");
	g_UMass_ampintegral_normalized->SetTitle("UMass");
	g_all_ampintegral_normalized->SetTitle("");

	double xmin = -25;
	double xmax = 25;
	double ymin = 0.8;
	double ymax = 2.4;

	g_INFN_blumlein->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_blumlein->GetXaxis()->SetLimits(xmin,xmax);
	g_all_blumlein->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_amp30->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_amp30->GetXaxis()->SetLimits(xmin,xmax);
	g_all_amp30->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_ampexp->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_ampexp->GetXaxis()->SetLimits(xmin,xmax);
	g_all_ampexp->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_ampintegral->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_ampintegral->GetXaxis()->SetLimits(xmin,xmax);
	g_all_ampintegral->GetXaxis()->SetLimits(xmin,xmax);

	g_INFN_blumlein_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_blumlein_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_all_blumlein_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_amp30_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_amp30_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_all_amp30_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_ampexp_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_ampexp_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_all_ampexp_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_INFN_ampintegral_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_UMass_ampintegral_normalized->GetXaxis()->SetLimits(xmin,xmax);
	g_all_ampintegral_normalized->GetXaxis()->SetLimits(xmin,xmax);

	g_INFN_blumlein_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_UMass_blumlein_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_all_blumlein_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_INFN_amp30_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_UMass_amp30_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_all_amp30_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_INFN_ampexp_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_UMass_ampexp_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_all_ampexp_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_INFN_ampintegral_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_UMass_ampintegral_normalized->GetYaxis()->SetRangeUser(ymin,ymax);
	g_all_ampintegral_normalized->GetYaxis()->SetRangeUser(ymin,ymax);


	g_INFN_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_UMass_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_all_blumlein->GetXaxis()->SetTitle("x [mm]");
	g_INFN_amp30->GetXaxis()->SetTitle("x [mm]");
	g_UMass_amp30->GetXaxis()->SetTitle("x [mm]");
	g_all_amp30->GetXaxis()->SetTitle("x [mm]");
	g_INFN_ampexp->GetXaxis()->SetTitle("x [mm]");
	g_UMass_ampexp->GetXaxis()->SetTitle("x [mm]");
	g_all_ampexp->GetXaxis()->SetTitle("x [mm]");
	g_INFN_ampintegral->GetXaxis()->SetTitle("x [mm]");
	g_UMass_ampintegral->GetXaxis()->SetTitle("x [mm]");
	g_all_ampintegral->GetXaxis()->SetTitle("x [mm]");

	g_INFN_blumlein_normalized->GetXaxis()->SetTitle("x [mm]");
	g_UMass_blumlein_normalized->GetXaxis()->SetTitle("x [mm]");
	g_all_blumlein_normalized->GetXaxis()->SetTitle("x [mm]");
	g_INFN_amp30_normalized->GetXaxis()->SetTitle("x [mm]");
	g_UMass_amp30_normalized->GetXaxis()->SetTitle("x [mm]");
	g_all_amp30_normalized->GetXaxis()->SetTitle("x [mm]");
	g_INFN_ampexp_normalized->GetXaxis()->SetTitle("x [mm]");
	g_UMass_ampexp_normalized->GetXaxis()->SetTitle("x [mm]");
	g_all_ampexp_normalized->GetXaxis()->SetTitle("x [mm]");
	g_INFN_ampintegral_normalized->GetXaxis()->SetTitle("x [mm]");
	g_UMass_ampintegral_normalized->GetXaxis()->SetTitle("x [mm]");
	g_all_ampintegral_normalized->GetXaxis()->SetTitle("x [mm]");


	g_INFN_blumlein->GetYaxis()->SetTitle("Blumlein [mG]");
	g_UMass_blumlein->GetYaxis()->SetTitle("Blumlein [mG]");
	g_all_blumlein->GetYaxis()->SetTitle("Blumlein [mG]");
	g_INFN_amp30->GetYaxis()->SetTitle("Transient at 30 #mus [mG]");
	g_UMass_amp30->GetYaxis()->SetTitle("Transient at 30 #mus [mG]");
	g_all_amp30->GetYaxis()->SetTitle("Transient at 30 #mus [mG]");
	g_INFN_ampexp->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [mG]");
	g_UMass_ampexp->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [mG]");
	g_all_ampexp->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [mG]");
	g_INFN_ampintegral->GetYaxis()->SetTitle("Transient integral [30,700] #mus [mG]");
	g_UMass_ampintegral->GetYaxis()->SetTitle("Transient integral [30,700] #mus [mG]");
	g_all_ampintegral->GetYaxis()->SetTitle("Transient integral [30,700] #mus [mG]");

	g_INFN_blumlein_normalized->GetYaxis()->SetTitle("Blumlein [arb. u.]");
	g_UMass_blumlein_normalized->GetYaxis()->SetTitle("Blumlein [arb. u.]");
	g_all_blumlein_normalized->GetYaxis()->SetTitle("Blumlein [arb. u.]");
	g_INFN_amp30_normalized->GetYaxis()->SetTitle("Transient at 30 #mus [arb. u.]");
	g_UMass_amp30_normalized->GetYaxis()->SetTitle("Transient at 30 #mus [arb. u.]");
	g_all_amp30_normalized->GetYaxis()->SetTitle("Transient at 30 #mus [arb. u.]");
	g_INFN_ampexp_normalized->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [arb. u.]");
	g_UMass_ampexp_normalized->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [arb. u.]");
	g_all_ampexp_normalized->GetYaxis()->SetTitle("Transient amplitude at 30 #mus [arb. u.]");
	g_INFN_ampintegral_normalized->GetYaxis()->SetTitle("Transient integral [30,700] #mus [arb. u.]");
	g_UMass_ampintegral_normalized->GetYaxis()->SetTitle("Transient integral [30,700] #mus [arb. u.]");
	g_all_ampintegral_normalized->GetYaxis()->SetTitle("Transient integral [30,700] #mus [arb. u.]");



	//Get umass model
	TFile* f_umass = TFile::Open("../Bk_calculation/UMass/UMass_model_0p1.root");
	TH2D* h2_umass = (TH2D*)f_umass->Get("h2");
	TH1D* h1_umass_y0 = (TH1D*)f_umass->Get("h1_y0");
	TH1D* h1_umass_y0_INFN = (TH1D*)h1_umass_y0->Clone("h1_umass_y0_INFN");
	TH1D* h1_umass_y0_UMass = (TH1D*)h1_umass_y0->Clone("h1_umass_y0_UMass");
	h1_umass_y0_INFN->Reset();
	h1_umass_y0_UMass->Reset();
	for (int bn=1; bn<=h1_umass_y0_INFN->GetXaxis()->GetNbins(); bn++){
		double biny_low = h2_umass->GetYaxis()->FindBin(-16);
		double biny_high = h2_umass->GetYaxis()->FindBin(16);
		double y = h2_umass->Integral(bn,bn,biny_low,biny_high);
		h1_umass_y0_INFN->SetBinContent(bn,y);
	}
	h1_umass_y0_INFN->Scale(1./h1_umass_y0_INFN->Interpolate(0));
	for (int bn=1; bn<=h1_umass_y0_UMass->GetXaxis()->GetNbins(); bn++){
		double biny_low = h2_umass->GetYaxis()->FindBin(-14.5);
		double biny_high = h2_umass->GetYaxis()->FindBin(14.5);
		double y = h2_umass->Integral(bn,bn,biny_low,biny_high);
		h1_umass_y0_UMass->SetBinContent(bn,y);
	}
	h1_umass_y0_UMass->Scale(1./h1_umass_y0_UMass->Interpolate(0));

	// now draw and do parabolic fits
	cout<<"\nDrawing \n";

	gROOT->SetBatch(kFALSE);
	gStyle->SetOptStat(0);


	TF1* f_quadratic = new TF1("f_quadratic","[0]+[1]*x*x");


	TGraph* g_normpoint = (TGraph*)g_INFN_ampexp_normalized->Clone("g_normpoint");
	g_normpoint->Set(0);
	g_normpoint->SetPoint(0,0,1);
	g_normpoint->SetMarkerStyle(20);
	g_normpoint->SetMarkerColor(7);

	TCanvas* can = new TCanvas("","",800,800);
	//can->Divide(2,1);
	//can->cd(1);
	//g_all_blumlein_normalized->Draw("APZ");
	//g_INFN_blumlein_normalized->Draw("PZ");
	//g_UMass_blumlein_normalized->Draw("PZ");
	//f_quadratic->SetParameters(1,0.001);
	//f_quadratic->SetLineColor(kViolet);
	//g_all_blumlein_normalized->Fit(f_quadratic);
	//gPad->SetGridx();
	//gPad->SetGridy();
	//TLegend* leg1 = new TLegend(0.35,0.6,0.65,0.8);
	//leg1->AddEntry(g_INFN_blumlein_normalized,"INFN","PL");
	//leg1->AddEntry(g_UMass_blumlein_normalized,"UMass","PL");
	//leg1->AddEntry(f_quadratic,Form("%.2f+%.5fx^{2}",f_quadratic->GetParameter(0),f_quadratic->GetParameter(1)),"L");
	//leg1->Draw();
	//g_normpoint->Draw("P");
	//can->cd(2);
	g_all_ampexp_normalized->Draw("APZ");
	g_INFN_ampexp_normalized->Draw("PZ");
	g_UMass_ampexp_normalized->Draw("PZ");
	h1_umass_y0->SetLineColor(kGreen+2);
	h1_umass_y0_INFN->SetLineColor(kGreen+2);
	h1_umass_y0_UMass->SetLineColor(kRed);
	h1_umass_y0->SetLineWidth(2);
	h1_umass_y0_INFN->SetLineWidth(2);
	h1_umass_y0_UMass->SetLineWidth(2);
	//h1_umass_y0->Draw("HIST L SAME");
	h1_umass_y0_INFN->Draw("HIST L SAME");
	//h1_umass_y0_UMass->Draw("HIST L SAME");
	//f_quadratic->SetParameters(1,0.005);
	//g_all_ampexp_normalized->Fit(f_quadratic);
	//f_quadratic->SetLineColor(kViolet);
	gPad->SetGridx();
	gPad->SetGridy();
	TLegend* leg2 = new TLegend(0.35,0.6,0.65,0.8);
	leg2->AddEntry(g_INFN_ampexp_normalized,"INFN","PL");
	leg2->AddEntry(g_UMass_ampexp_normalized,"UMass","PL");
	//leg2->AddEntry(f_quadratic,Form("%.2f+%.5fx^{2}",f_quadratic->GetParameter(0),f_quadratic->GetParameter(1)),"L");
	leg2->AddEntry(h1_umass_y0_INFN,"UMass radial model","L");
	//leg2->AddEntry(h1_umass_y0_INFN,"UMass radial model","L");
	//leg2->AddEntry(h1_umass_y0_UMass,"UMass radial model","L");
	leg2->Draw();
	g_normpoint->Draw("P");


	TFile* fout = new TFile("paper_radial_plot.root","recreate");
	g_all_ampexp_normalized->Write("g_all");
	g_INFN_ampexp_normalized->Write("g_INFN");
	g_UMass_ampexp_normalized->Write("g_UMass");
	h1_umass_y0_INFN->Write("model");
	g_normpoint->Write("normalization_point");
	fout->Write();
	fout->Close();

}
#include "../analysis_tools.C"

void create_INFN_root(){

	//R0

	TFile* f1_p = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f1_n = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan28_B5173_H25_Q00.root");
	TFile* f1_0 = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");
	TH1D* h1_kick1_R0_p = ((TProfile*)f1_p->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_n = ((TProfile*)f1_n->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_0 = ((TProfile*)f1_0->Get("trace_kick1"))->ProjectionX();


	double blum_norm_x = -0.323;
	double blum_norm_y_R0 = 129.;
	double blum_norm_y_R1 = 157.;
	double blum_norm_x_trace = 4.769;
	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1_kick1_R0 = (TH1D*)h1_kick1_R0_p->Clone("h1_kick1_R0");
	h1_kick1_R0->Scale(0.5);
	h1_kick1_R0->Add(h1_kick1_R0_n,-0.5);
	h1_kick1_R0->Add(h1_kick1_R0_0,0.2);
	cleanTrace(h1_kick1_R0,-200);
	TH1D* h1_kick1_R0_ra = smoothing(h1_kick1_R0,"");

	cleanTrace(h1_kick1_R0_p,-200);
	TH1D* h1_kick1_R0_p_ra = smoothing(h1_kick1_R0_p,"");

	double r0_norm = blum_norm_y_R0/h1_kick1_R0_ra->Interpolate(blum_norm_x);
	h1_kick1_R0->Scale(r0_norm);
	h1_kick1_R0_ra->Scale(r0_norm);

	r0_norm = blum_norm_y_R0/h1_kick1_R0_p_ra->Interpolate(blum_norm_x);
	h1_kick1_R0_p->Scale(r0_norm);
	h1_kick1_R0_p_ra->Scale(r0_norm);



	//R1

	TFile* f_R1 = TFile::Open("../Eddy_analysis/analysis/analysis_SD_R1_eddy_oct8_H0_nofilter_Bfield.root");
	TH1D* h1_kick1_R1 = ((TProfile*)f_R1->Get("trace_kick1"))->ProjectionX();


	cleanTrace(h1_kick1_R1,-200);
	TH1D* h1_kick1_R1_ra = smoothing(h1_kick1_R1,"");

	double r1_norm = blum_norm_y_R1/h1_kick1_R1_ra->Interpolate(blum_norm_x);
	h1_kick1_R1_ra->Scale(r1_norm);



	TFile* fout = new TFile("INFN.root","recreate");


	h1_kick1_R0_p->SetTitle("INFN K3 0.0 mm");
	h1_kick1_R0_p_ra->SetTitle("INFN K3 0.0 mm");
	h1_kick1_R0->SetTitle("INFN K3 0.0 mm");
	h1_kick1_R0_ra->SetTitle("INFN K3 0.0 mm");
	h1_kick1_R1->SetTitle("INFN K3 17.5 mm");
	h1_kick1_R1_ra->SetTitle("INFN K3 17.5 mm");

	h1_kick1_R0_p->Write("h1_kick1_R0_p");
	h1_kick1_R0_p_ra->Write("h1_kick1_R0_p_ra");
	h1_kick1_R0->Write("h1_kick1_R0");
	h1_kick1_R0_ra->Write("h1_kick1_R0_ra");
	h1_kick1_R1->Write("h1_kick1_R1");
	h1_kick1_R1_ra->Write("h1_kick1_R1_ra");


	fout->Close();



}
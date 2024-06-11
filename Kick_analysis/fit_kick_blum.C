#include "../analysis_tools.C"

TGraph* h_kick_template;


double template_function(double* x, double* pars){
	double xx = x[0];
	double A = pars[0];
	double t0 = pars[1];

	double t_graph = 1.e3*(xx-t0);
	if (t_graph < -2 || t_graph > 6) return 0;
	return A*h_kick_template->Eval(t_graph);
}


void fit_kick_blum(){

	TFile* f0 = TFile::Open("output_FD_R0_oct9_H0_Bfield.root");
	TFile* f1 = TFile::Open("output_FD_R1_oct8_H0_Bfield.root");

	TFile* f0_blum = TFile::Open("analysis_FD_R0_blum_oct9_H22p5_Bfield.root");
	TFile* f1_blum = TFile::Open("analysis_FD_R1_blum_oct8_H0_Bfield.root");

	TGraph*	h1_kick1_R0_norm = (TGraph*)f0->Get("normalized_kick_1");
	TGraph*	h1_kick1_R1_norm = (TGraph*)f1->Get("normalized_kick_1");
	h1_kick1_R0_norm->SetTitle("Kick R0 (normalized)");
	h1_kick1_R1_norm->SetTitle("Kick R1 (normalized)");

	TH1D* h1_kick1_R0_blum = ((TProfile*)f0_blum->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R1_blum = ((TProfile*)f1_blum->Get("trace_kick1"))->ProjectionX();
	h1_kick1_R1_blum->Scale(-1);

	TH1D* h1_kick1_R0_blum_zoom = (TH1D*)h1_kick1_R0_blum->Clone("h1_kick1_R0_blum_zoom");
	TH1D* h1_kick1_R1_blum_zoom = (TH1D*)h1_kick1_R1_blum->Clone("h1_kick1_R1_blum_zoom");

	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;
	TH1D* h1_kick1_R0_blum_ra = runningAverage(runningAverage(runningAverage(h1_kick1_R0_blum,100,true),200,true),300,true,"");
	TH1D* h1_kick1_R1_blum_ra = runningAverage(runningAverage(runningAverage(h1_kick1_R1_blum,100,true),200,true),300,true,"");


	h_kick_template = h1_kick1_R0_norm;

	TF1* f_template = new TF1("f_template",template_function,-0.001,0.01,2);
	f_template->SetNpx(1000);
	f_template->SetParNames("A","t0");

	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.45,-0.15);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);


	new TCanvas();
	h1_kick1_R0_blum->SetLineWidth(2);
	h1_kick1_R0_blum->GetXaxis()->SetRangeUser(-0.001,0.005);
	h1_kick1_R0_blum->GetYaxis()->SetRangeUser(-150,150);
	h1_kick1_R0_blum->Draw("HIST");

	f_template->SetParameters(1000,0);
	f_template->SetParLimits(0,1,1e7);
	f_template->SetParLimits(1,-5e-4,5e-4);
	TFitResultPtr fit_R0_kick = h1_kick1_R0_blum->Fit(f_template,"S","",0.0014,0.005);
	f_template->DrawCopy("SAME");

	double fit_kick_R0 = fit_R0_kick->Parameter(0);
	double fit_kick_R0_err = fit_R0_kick->ParError(0);

	new TCanvas();
	h1_kick1_R1_blum->SetLineWidth(2);
	h1_kick1_R1_blum->GetXaxis()->SetRangeUser(-0.001,0.005);
	h1_kick1_R1_blum->GetYaxis()->SetRangeUser(-150,150);
	h1_kick1_R1_blum->Draw("HIST");

	f_template->SetParameters(2000,0);
	f_template->SetParLimits(0,1,1e7);
	f_template->SetParLimits(1,-5e-4,5e-4);
	TFitResultPtr fit_R1_kick = h1_kick1_R1_blum->Fit(f_template,"S","",0.001,0.005);
	f_template->DrawCopy("SAME");

	double fit_kick_R1 = fit_R1_kick->Parameter(0);
	double fit_kick_R1_err = fit_R1_kick->ParError(0);



	f_blumlein->SetParameters(0.5,-0.3,-1000.0);
	new TCanvas();
	h1_kick1_R0_blum_zoom->GetXaxis()->SetRangeUser(-0.6,0.1);
	h1_kick1_R0_blum_zoom->GetYaxis()->SetRangeUser(-2,2);
	h1_kick1_R0_blum_zoom->DrawCopy("HIST");
	TFitResultPtr fit_R0_blum = h1_kick1_R0_blum_zoom->Fit("f_blumlein","S","",-0.45,-0.15);//-0.4,-0.2);
	f_blumlein->DrawCopy("SAME");

	double fit_blum_R0 = abs(fit_R0_blum->Parameter(0));
	double fit_blum_R0_err = fit_R0_blum->ParError(0);

	double ratio_R0 = fit_kick_R0/fit_blum_R0;
	double ratio_R0_err = ratio_R0*sqrt((fit_blum_R0_err/fit_blum_R0)*(fit_blum_R0_err/fit_blum_R0) + (fit_kick_R0_err/fit_kick_R0)*(fit_kick_R0_err/fit_kick_R0));


	f_blumlein->SetParameters(0.5,-0.3,-1000.0);
	new TCanvas();
	h1_kick1_R1_blum_zoom->GetXaxis()->SetRangeUser(-0.6,0.1);
	h1_kick1_R1_blum_zoom->GetYaxis()->SetRangeUser(-2,2);
	h1_kick1_R1_blum_zoom->DrawCopy("HIST");
	TFitResultPtr fit_R1_blum = h1_kick1_R1_blum_zoom->Fit("f_blumlein","S","",-0.45,-0.15);//-0.4,-0.2);
	f_blumlein->DrawCopy("SAME");

	double fit_blum_R1 = abs(fit_R1_blum->Parameter(0));
	double fit_blum_R1_err = fit_R1_blum->ParError(0);

	double ratio_R1 = fit_kick_R1/fit_blum_R1;
	double ratio_R1_err = ratio_R1*sqrt((fit_blum_R1_err/fit_blum_R1)*(fit_blum_R1_err/fit_blum_R1) + (fit_kick_R1_err/fit_kick_R1)*(fit_kick_R1_err/fit_kick_R1));


	cout<<"R0 ---\n";
	cout<<"blum: "<<fit_blum_R0<<" +- "<<fit_blum_R0_err<<" mV\n";
	cout<<"kick: "<<fit_kick_R0<<" +- "<<fit_kick_R0_err<<" mV\n";
	cout<<"Ratio: "<<ratio_R0<<" +- "<<ratio_R0_err<<"\n";

	cout<<"R1 ---\n";
	cout<<"blum: "<<fit_blum_R1<<" +- "<<fit_blum_R1_err<<" mV\n";
	cout<<"kick: "<<fit_kick_R1<<" +- "<<fit_kick_R1_err<<" mV\n";
	cout<<"Ratio: "<<ratio_R1<<" +- "<<ratio_R1_err<<"\n";

}
#include "../analysis_tools.C"

void fit_golden_trace(){


	TFile* f1_p = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f1_n = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan28_B5173_H25_Q00.root");
	TFile* f1_0 = TFile::Open("../Eddy_analysis/analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");
	TH1D* h1_kick1_R0_p = ((TProfile*)f1_p->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_n = ((TProfile*)f1_n->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_R0_0 = ((TProfile*)f1_0->Get("trace_kick1"))->ProjectionX();

	cleanTrace(h1_kick1_R0_p,-500);
	cleanTrace(h1_kick1_R0_n,-500);
	cleanTrace(h1_kick1_R0_0,-500);
	cleanTrace(h1_kick1_R0_p,500);
	cleanTrace(h1_kick1_R0_n,500);
	cleanTrace(h1_kick1_R0_0,500);

	double blum_norm_x = -0.323;
	double blum_norm_y_R0 = 129.;
	double blum_norm_y_R1 = 157.;
	double blum_norm_x_trace = 4.769;
	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;


	TH1D* h1_kick1_R0 = (TH1D*)h1_kick1_R0_p->Clone("h1_kick1_R0");
	h1_kick1_R0->Scale(0.5);
	h1_kick1_R0->Add(h1_kick1_R0_n,-0.5);
	h1_kick1_R0->Add(h1_kick1_R0_0,0.2);
	h1_kick1_R0->GetYaxis()->SetTitle("B field [mG]");
	cleanTrace(h1_kick1_R0,-60);

	TH1D* h1_kick1_R0_ra = smoothing(h1_kick1_R0,"");
	double r0_norm = blum_norm_y_R0/h1_kick1_R0_ra->Interpolate(blum_norm_x);
	h1_kick1_R0_ra->Scale(r0_norm);
	h1_kick1_R0->Scale(r0_norm);

	TH1D* h1_kick1_R0_ra2 = runningAverage(h1_kick1_R0_ra,150,true,"");
	TH1D* h1_kick1_R0_ra_detrended = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_detrended");
	h1_kick1_R0_ra_detrended->Add(h1_kick1_R0_ra2,-1);

	TH1D* h1_FFT_raw = doFFT(h1_kick1_R0,0.5,1.0);
	TH1D* h1_FFT = doFFT(h1_kick1_R0_ra,0.5,1.0);


	double xmin = -0.1;
	double xmax = 1.0;
	double ymin = -20;
	double ymax = 10;

	double xmin_fit = 0.03;
	double xmax_fit = 0.7;




	TLine* l30 = new TLine(0.03,ymin,0.03,ymax);
	TLine* l700 = new TLine(0.7,ymin,0.7,ymax);
	TLine* l0h = new TLine(xmin,0,xmax,0);
	TLine* l0v = new TLine(0,ymin,0,ymax);
	l30->SetLineStyle(kSolid);
	l30->SetLineColor(kGreen);
	l30->SetLineWidth(2);
	l700->SetLineStyle(kSolid);
	l700->SetLineColor(kGreen);
	l700->SetLineWidth(2);
	l0h->SetLineStyle(kDashed);
	l0h->SetLineWidth(2);
	l0v->SetLineStyle(kDashed);
	l0v->SetLineWidth(2);

	TGraph* g_box_30 = new TGraph(5);
	g_box_30->SetPoint(0,xmin,ymin);
	g_box_30->SetPoint(1,xmin_fit,ymin);
	g_box_30->SetPoint(2,xmin_fit,ymax);
	g_box_30->SetPoint(3,xmin,ymax);
	g_box_30->SetPoint(4,xmin,ymin);
	g_box_30->SetFillStyle(3144);

	TGraph* g_box_700 = new TGraph(5);
	g_box_700->SetPoint(0,xmax,ymin);
	g_box_700->SetPoint(1,xmax_fit,ymin);
	g_box_700->SetPoint(2,xmax_fit,ymax);
	g_box_700->SetPoint(3,xmax,ymax);
	g_box_700->SetPoint(4,xmax,ymin);
	g_box_700->SetFillStyle(3144);





	//Drawing

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);

	new TCanvas();
	h1_kick1_R0_ra->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_kick1_R0_ra->SetLineWidth(2);
	h1_kick1_R0_ra->Draw("HIST");
	h1_kick1_R0_ra2->Draw("HIST SAME");
	h1_kick1_R0_ra_detrended->Draw("HIST SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");


	new TCanvas();
	h1_kick1_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_kick1_R0->SetLineColor(kRed);
	h1_kick1_R0->Draw("HIST");
	h1_kick1_R0_ra->Draw("HIST SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");

	new TCanvas();
	h1_FFT_raw->SetLineColor(kRed);
	h1_FFT->Draw("HIST");
	h1_FFT_raw->Draw("HIST SAME");


	//Fitting

	TF1* f1 = new TF1("f1","[0]-[1]*exp(-x/[2])",xmin_fit,xmax_fit);
	f1->SetParameters(0,10,0.07);
	f1->SetParLimits(1,5,100);
	f1->SetParLimits(2,0.02,0.2);
	f1->SetParNames("b","A","#tau");
	f1->FixParameter(0,0.0);

	TF1* f1bis = new TF1("f1bis","[0]-[1]*exp(-x/[2])*sin([3]*6.283185307*x+[4])",xmin_fit,xmax_fit);
	f1bis->SetParameters(0,10,0.1,2,0);
	f1bis->SetParLimits(1,5,100);
	f1bis->SetParLimits(2,0.01,0.2);
	f1bis->SetParLimits(3,1,5);
	f1bis->FixParameter(0,0.0);
	f1bis->SetParNames("b","A","#tau","f","#phi");

	TF1* f1bis_exp = new TF1("f1bis_exp","[0]-[1]*exp(-x/[2])",xmin_fit,xmax_fit);
	f1bis_exp->SetParNames("b","A","#tau");

	TF1* f2 = new TF1("f2","[0]-[1]*exp(-x/[2])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])",xmin_fit,xmax_fit);
	f2->SetParameters(0,10,0.07,2,0.3,17,0);
	f2->SetParLimits(1,5,100);
	f2->SetParLimits(2,0.02,0.2);
	f2->SetParLimits(3,0.5,10);
	f2->SetParLimits(4,0.1,1.0);
	f2->SetParLimits(5,14,19);
	f2->FixParameter(0,0.0);
	f2->SetParNames("b","A","#tau","A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	TF1* f2_exp = new TF1("f2_exp","[0]-[1]*exp(-x/[2])",xmin_fit,xmax_fit);
	f2_exp->SetParNames("b","A","#tau");
	TF1* f2_sine = new TF1("f2_sine","[0]*exp(-x/[1])*sin([2]*6.283185307*x+[3])",-0.1,xmax_fit);
	f2_sine->SetParNames("A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	TF1* f3 = new TF1("f3","[0]-[1]*exp(-x/[2])*sin([3]*6.283185307*x+[4])+[5]*exp(-x/[6])*sin([7]*6.283185307*x+[8])",xmin_fit,xmax_fit);
	f3->SetParameters(0,10,0.1,2,0,2,0.3,17,0);
	f3->SetParLimits(1,5,100);
	f3->SetParLimits(2,0.01,0.2);
	f3->SetParLimits(3,1,5);
	f3->SetParLimits(5,0.5,10);
	f3->SetParLimits(6,0.1,1.0);
	f3->SetParLimits(7,14,19);
	f3->FixParameter(0,0.0);
	f3->SetParNames("b","A","#tau","f","#phi","A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	TF1* f3_exp = new TF1("f3_exp","[0]-[1]*exp(-x/[2])*sin([3]*6.283185307*x+[4])",xmin_fit,xmax_fit);
	f3_exp->SetParNames("b","A","#tau","f","#phi");
	TF1* f3_sine = new TF1("f3_sine","[0]*exp(-x/[1])*sin([2]*6.283185307*x+[3])",-0.1,xmax_fit);
	f3_sine->SetParNames("A_{1}","#tau_{1}","f_{1}","#phi_{1}");


	TF1* f4 = new TF1("f4","[0]-[1]*exp(-x/[2])*sin([3]*6.283185307*x+[4])+[5]*exp(-x/[6])*sin([7]*6.283185307*x+[8])+[9]*sin([10]*6.283185307*x+[11])",xmin_fit,xmax_fit);
	f4->SetParameters(0,10,0.1,2,0,2,0.3,17,0);
	f4->SetParameters(9,1);
	f4->SetParameters(10,6);
	f4->SetParameters(11,0);
	f4->SetParLimits(1,5,100);
	f4->SetParLimits(2,0.01,0.2);
	f4->SetParLimits(3,1,5);
	f4->SetParLimits(5,0.5,10);
	f4->SetParLimits(6,0.1,1.0);
	f4->SetParLimits(7,14,19);
	f4->SetParLimits(8,0.1,10);
	f4->SetParLimits(9,0.1,5);
	f4->SetParLimits(10,4,10);
	f4->FixParameter(0,0.0);
	f4->SetParNames("b","A","#tau","f","#phi","A_{1}","#tau_{1}","f_{1}","#phi_{1}","A_{2}","f_{2}");
	f4->SetParName(11,"#phi_{2}");

	TF1* f4_exp = new TF1("f4_exp","[0]-[1]*exp(-x/[2])*sin([3]*6.283185307*x+[4])",xmin_fit,xmax_fit);
	f4_exp->SetParNames("b","A","#tau","f","#phi");
	TF1* f4_sine = new TF1("f4_sine","[0]*exp(-x/[1])*sin([2]*6.283185307*x+[3])",-0.1,xmax_fit);
	f4_sine->SetParNames("A_{1}","#tau_{1}","f_{1}","#phi_{1}");
	TF1* f4_sine2 = new TF1("f4_sine2","[0]*sin([1]*6.283185307*x+[2])",xmin_fit,xmax_fit);
	f4_sine2->SetParNames("A_{2}","f_{2}","#phi_{2}");
	

	TFitResultPtr fit_res1 = h1_kick1_R0_ra->Fit(f1,"SN0","",0.03,0.7);
	TH1D* h1_kick1_R0_ra_res1 = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_res1");
	h1_kick1_R0_ra_res1->SetTitle("Fit residual");
	h1_kick1_R0_ra_res1->Add(f1,-1);
	h1_kick1_R0_ra_res1->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra_res1->GetYaxis()->SetRangeUser(ymin,ymax);


	f1bis->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
	TFitResultPtr fit_res1bis = h1_kick1_R0_ra->Fit(f1bis,"SN0","",0.03,0.7);
	TH1D* h1_kick1_R0_ra_res1bis = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_res1bis");
	h1_kick1_R0_ra_res1bis->SetTitle("Fit residual");
	h1_kick1_R0_ra_res1bis->Add(f1bis,-1);
	h1_kick1_R0_ra_res1bis->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra_res1bis->GetYaxis()->SetRangeUser(ymin,ymax);


	f2->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
	TFitResultPtr fit_res2 = h1_kick1_R0_ra->Fit(f2,"SN0","",0.03,0.7);
	f2_exp->SetParameters(f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2));
	f2_sine->SetParameters(f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5),f2->GetParameter(6));
	TH1D* h1_kick1_R0_ra_res2 = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_res2");
	h1_kick1_R0_ra_res2->SetTitle("Fit residual");
	h1_kick1_R0_ra_res2->Add(f2,-1);
	h1_kick1_R0_ra_res2->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra_res2->GetYaxis()->SetRangeUser(ymin,ymax);
	TH1D* h1_fft_res2 = doFFT(h1_kick1_R0_ra_res2,xmin_fit,xmax_fit);


	f3->SetParameters(f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5),f2->GetParameter(6));
	TFitResultPtr fit_res3 = h1_kick1_R0_ra->Fit(f3,"SN0","",0.03,0.7);
	f3_exp->SetParameters(f3->GetParameter(0),f3->GetParameter(1),f3->GetParameter(2),f3->GetParameter(3),f3->GetParameter(4));
	f3_sine->SetParameters(f3->GetParameter(5),f3->GetParameter(6),f3->GetParameter(7),f3->GetParameter(8));
	TH1D* h1_kick1_R0_ra_res3 = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_res3");
	h1_kick1_R0_ra_res3->SetTitle("Fit residual");
	h1_kick1_R0_ra_res3->Add(f3,-1);
	h1_kick1_R0_ra_res3->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra_res3->GetYaxis()->SetRangeUser(ymin,ymax);
	TH1D* h1_fft_res3 = doFFT(h1_kick1_R0_ra_res3,xmin_fit,xmax_fit);


	f4->SetParameters(f3->GetParameter(0),f3->GetParameter(1),f3->GetParameter(2),f3->GetParameter(3),f3->GetParameter(4),f3->GetParameter(5),f3->GetParameter(6),f3->GetParameter(7),f3->GetParameter(8));
	TFitResultPtr fit_res4 = h1_kick1_R0_ra->Fit(f4,"SN0","",0.03,0.7);
	f4_exp->SetParameters(f4->GetParameter(0),f4->GetParameter(1),f4->GetParameter(2),f4->GetParameter(3),f4->GetParameter(4));
	f4_sine->SetParameters(f4->GetParameter(5),f4->GetParameter(6),f4->GetParameter(7),f4->GetParameter(8));
	f4_sine2->SetParameters(f4->GetParameter(9),f4->GetParameter(10),f4->GetParameter(11));
	TH1D* h1_kick1_R0_ra_res4 = (TH1D*)h1_kick1_R0_ra->Clone("h1_kick1_R0_ra_res4");
	h1_kick1_R0_ra_res4->SetTitle("Fit residual");
	h1_kick1_R0_ra_res4->Add(f4,-1);
	h1_kick1_R0_ra_res4->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_kick1_R0_ra_res4->GetYaxis()->SetRangeUser(ymin,ymax);
	TH1D* h1_fft_res4 = doFFT(h1_kick1_R0_ra_res4,xmin_fit,xmax_fit);


	TCanvas* can1 = new TCanvas("can1","",1800,800);
	can1->Divide(2,1);
	can1->cd(1);
	h1_kick1_R0_ra->Draw("HIST");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");
	f1->Draw("SAME");
	can1->cd(2);
	h1_kick1_R0_ra_res1->Draw("HIST");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");


	TCanvas* can1bis = new TCanvas("can1bis","",1800,800);
	can1bis->Divide(2,1);
	can1bis->cd(1);
	h1_kick1_R0_ra->Draw("HIST");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");
	f1bis->Draw("SAME");
	can1bis->cd(2);
	h1_kick1_R0_ra_res1bis->Draw("HIST");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");


	TCanvas* can2 = new TCanvas("can2","",1800,800);
	can2->Divide(2,1);
	can2->cd(1);
	h1_kick1_R0_ra->Draw("HIST");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");
	f2->Draw("SAME");
	can2->cd(2);
	h1_kick1_R0_ra_res2->Draw("HIST");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");
	new TCanvas();
	h1_fft_res2->Draw("HIST");
	new TCanvas();
	h1_kick1_R0_ra->Draw("HIST");
	f2_exp->Draw("SAME");
	f2_sine->Draw("SAME");

	TCanvas* can3 = new TCanvas("can3","",1800,800);
	can3->Divide(2,1);
	can3->cd(1);
	h1_kick1_R0_ra->Draw("HIST");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");
	f3->Draw("SAME");
	can3->cd(2);
	h1_kick1_R0_ra_res3->Draw("HIST");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");
	new TCanvas();
	h1_fft_res3->Draw("HIST");
	new TCanvas();
	h1_kick1_R0_ra->Draw("HIST");
	f3_exp->Draw("SAME");
	f3_sine->Draw("SAME");


	TCanvas* can4 = new TCanvas("can4","",1800,800);
	can4->Divide(2,1);
	can4->cd(1);
	h1_kick1_R0_ra->Draw("HIST");
	l30->Draw("SAME");
	l700->Draw("SAME");
	l0h->Draw("SAME");
	l0v->Draw("SAME");
	f4->Draw("SAME");
	can4->cd(2);
	h1_kick1_R0_ra_res4->Draw("HIST");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");
	new TCanvas();
	h1_fft_res4->Draw("HIST");
	new TCanvas();
	h1_kick1_R0_ra->Draw("HIST");
	f4_exp->Draw("SAME");
	f4_sine->Draw("SAME");
	f4_sine2->Draw("SAME");

	new TCanvas();
	h1_kick1_R0_ra_res1->Draw("HIST");
	h1_kick1_R0_ra_res1bis->Draw("HIST");
	h1_kick1_R0_ra_res2->Draw("HIST SAME");
	h1_kick1_R0_ra_res3->Draw("HIST SAME");
	h1_kick1_R0_ra_res4->Draw("HIST SAME");
	l0h->Draw("SAME");
	l30->Draw("SAME");
	l700->Draw("SAME");
	g_box_30->Draw("F");
	g_box_700->Draw("F");

	can1->SaveAs("fit_f1.png");
	can1bis->SaveAs("fit_f1bis.png");
	//can2->SaveAs("fit_f2.png");
	can3->SaveAs("fit_f3.png");
	can4->SaveAs("fit_f4.png");

	cout<<fit_res1->Chi2()<<" "<<fit_res1->Ndf()<<"\n";
	cout<<fit_res1bis->Chi2()<<" "<<fit_res1bis->Ndf()<<"\n";
	cout<<fit_res2->Chi2()<<" "<<fit_res2->Ndf()<<"\n";
	cout<<fit_res3->Chi2()<<" "<<fit_res3->Ndf()<<"\n";
	cout<<fit_res4->Chi2()<<" "<<fit_res4->Ndf()<<"\n";


	cout<<f1->GetExpFormula()<<"\n";
	cout<<f1bis->GetExpFormula()<<"\n";
	cout<<f3->GetExpFormula()<<"\n";
	cout<<f4->GetExpFormula()<<"\n";
}
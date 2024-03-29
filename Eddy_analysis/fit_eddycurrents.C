#include "../analysis_tools.C"

void fit_eddycurrents(TString input_file, TString output_file="", int useSNR=2){

	TFile* f = TFile::Open(input_file);

	TFile* fout;

	bool saveOutput = false;
	if (output_file != ""){
		saveOutput = true;
		fout = new TFile(output_file,"recreate");
	}

	TString SNRstring = useSNR < 0 ? "" : Form("_SNR%d",useSNR);
	cout<<"Using SNR string: "<<SNRstring<<"\n";

	TH1D* trace = ((TProfile*)f->Get(Form("trace%s",SNRstring.Data())))->ProjectionX();
	TH1D** kicks = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks[i] = ((TProfile*)f->Get(Form("trace%s_kick%d_calibrated",SNRstring.Data(),i+1)))->ProjectionX();
	}
	TH1D* kick8long = ((TProfile*)f->Get(Form("trace%s_kick8long",SNRstring.Data())))->ProjectionX();

	cleanTrace(trace,-50);
	for (int i=0; i<8; i++){
		cleanTrace(kicks[i],-50);
	}
	cleanTrace(kick8long,-50);

	TH1D* trace_ra = runningAverage_5_10_15(trace);
	trace_ra->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),51015));	
	TH1D** kicks_ra = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra[i] = runningAverage_5_10_15(kicks[i]);
		kicks_ra[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),51015));	
	}
	TH1D* kick8long_ra = runningAverage_5_10_15(kick8long);
	kick8long_ra->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),51015));	


	new TCanvas();
	trace->SetLineWidth(2);
	trace_ra->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_ra->SetLineColor(kRed);
	trace->Draw("HIST");
	trace_ra->Draw("HIST SAME");


	for (int i=0; i<8; i++){
		kicks[i]->GetYaxis()->SetRangeUser(-60,60);
		kicks_ra[i]->GetYaxis()->SetRangeUser(-60,60);
	}


	TCanvas* can_kicks = new TCanvas("can_kicks","",1600,800);
	can_kicks->Divide(2,1);
	can_kicks->cd(1);
	kicks[0]->SetLineWidth(2);
	kicks_ra[0]->SetLineWidth(2);
	kicks[0]->SetLineColor(kBlue);
	kicks_ra[0]->SetLineColor(kRed);
	kicks[0]->Draw("HIST");
	kicks_ra[0]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);
	can_kicks->cd(2);
	kicks[7]->SetLineWidth(2);
	kicks_ra[7]->SetLineWidth(2);
	kicks[7]->SetLineColor(kBlue);
	kicks_ra[7]->SetLineColor(kRed);
	kicks[7]->Draw("HIST");
	kicks_ra[7]->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.8,0.9,0.9);



	TH1D** h1_fft_kicks = new TH1D* [8];
	TH1D** h1_fft_kicks_ra = new TH1D* [8];
	double fft_xmin = 0.1; //ms from kick
	double fft_xmax = 8.0; 
	double dt = kicks[0]->GetBinWidth(1);
	for (int i=0; i<8; i++){
		h1_fft_kicks[i] = doFFT(kicks[i],fft_xmin,fft_xmax);
		h1_fft_kicks[i]->GetYaxis()->SetRangeUser(0,0.01);
	}
	for (int i=0; i<8; i++){
		h1_fft_kicks_ra[i] = doFFT(kicks_ra[i],fft_xmin,fft_xmax);
		h1_fft_kicks_ra[i]->GetYaxis()->SetRangeUser(0,0.01);
	}


	TText* text = new TText();
	TCanvas* can = new TCanvas("can","",1400,800);
	can->Divide(2,2);
	can->cd(1);
	kicks[0]->Draw("HIST");
	can->cd(2);
	kicks_ra[0]->Draw("HIST");
	can->cd(3);
	h1_fft_kicks[0]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks[0]->Integral()));
	gPad->SetLogx();
	can->cd(4);
	h1_fft_kicks_ra[0]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks_ra[0]->Integral()));
	gPad->SetLogx();


	TCanvas* can8 = new TCanvas("can8","",1400,800);
	can8->Divide(2,2);
	can8->cd(1);
	kicks[7]->Draw("HIST");
	can8->cd(2);
	kicks_ra[7]->Draw("HIST");
	can8->cd(3);
	h1_fft_kicks[7]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks[7]->Integral()));
	gPad->SetLogx();
	can8->cd(4);
	h1_fft_kicks_ra[7]->Draw("HIST");
	text->DrawTextNDC(0.6,0.8,Form("Integral: %.5f",0.5*h1_fft_kicks_ra[7]->Integral()));
	gPad->SetLogx();






	//Fits
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])");
	f_blumlein->SetParameters(40.0,-0.3,-1000.0);
	TFitResultPtr fit_blumlein = kicks_ra[0]->Fit("f_blumlein","QS+","",-0.4,-0.2);
	TFitResultPtr fit_baseline = kicks_ra[0]->Fit("pol0","QS+","",-1.0,-0.6);
	double peak=0;
	double peakerr=0;
	double base=0;
	double baseerr=0;
	if (fit_blumlein>=0){
		peak = fit_blumlein->Parameter(0);
		peakerr = fit_blumlein->ParError(0);
	}
	if (fit_baseline>=0){
		base = fit_baseline->Parameter(0);
		baseerr = fit_baseline->ParError(0);
	}
	double blumlein = peak - base;
	double blumerr = sqrt(peakerr*peakerr + baseerr*baseerr);
		
	cout<<"First blumlein: "<<blumlein<<" +- "<<blumerr<<"\n";




	TH1D** kicks_fit_exp = new TH1D* [8];
	TH1D** kicks_fit_exp_res = new TH1D* [8];
	TH1D** kicks_fit_exp17 = new TH1D* [8];
	TH1D** kicks_fit_exp17_res = new TH1D* [8];
	TH1D** kicks_fit_exp17_4 = new TH1D* [8];
	TH1D** kicks_fit_exp17_4_res = new TH1D* [8];
	TH1D** kicks_fit_exp17_4_20 = new TH1D* [8];
	TH1D** kicks_fit_exp17_4_20_res = new TH1D* [8];
	for (int i=0; i<8; i++){
		kicks_fit_exp[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp",i+1));
		kicks_fit_exp_res[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp_res",i+1));
		kicks_fit_exp17[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17",i+1));
		kicks_fit_exp17_res[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17_res",i+1));
		kicks_fit_exp17_4[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17_4",i+1));
		kicks_fit_exp17_4_res[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17_4_res",i+1));
		kicks_fit_exp17_4_20[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17_4_20",i+1));
		kicks_fit_exp17_4_20_res[i] = (TH1D*)kicks_ra[i]->Clone(Form("kick%d_fit_exp17_4_20_res",i+1));
		
		kicks_fit_exp[i]->SetLineColor(kBlue);
		kicks_fit_exp_res[i]->SetLineColor(kBlue);
		kicks_fit_exp[i]->SetLineWidth(1);
		kicks_fit_exp_res[i]->SetLineWidth(1);
		kicks_fit_exp[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp_res[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp[i]->GetYaxis()->SetRangeUser(-20,20);
		kicks_fit_exp_res[i]->GetYaxis()->SetRangeUser(-20,20);

		kicks_fit_exp17[i]->SetLineColor(kBlue);
		kicks_fit_exp17_res[i]->SetLineColor(kBlue);
		kicks_fit_exp17[i]->SetLineWidth(1);
		kicks_fit_exp17_res[i]->SetLineWidth(1);
		kicks_fit_exp17[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17_res[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17[i]->GetYaxis()->SetRangeUser(-20,20);
		kicks_fit_exp17_res[i]->GetYaxis()->SetRangeUser(-20,20);

		kicks_fit_exp17_4[i]->SetLineColor(kBlue);
		kicks_fit_exp17_4_res[i]->SetLineColor(kBlue);
		kicks_fit_exp17_4[i]->SetLineWidth(1);
		kicks_fit_exp17_4_res[i]->SetLineWidth(1);
		kicks_fit_exp17_4[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17_4_res[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17_4[i]->GetYaxis()->SetRangeUser(-20,20);
		kicks_fit_exp17_4_res[i]->GetYaxis()->SetRangeUser(-20,20);

		kicks_fit_exp17_4_20[i]->SetLineColor(kBlue);
		kicks_fit_exp17_4_20_res[i]->SetLineColor(kBlue);
		kicks_fit_exp17_4_20[i]->SetLineWidth(1);
		kicks_fit_exp17_4_20_res[i]->SetLineWidth(1);
		kicks_fit_exp17_4_20[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17_4_20_res[i]->GetXaxis()->SetRangeUser(0,0.8);
		kicks_fit_exp17_4_20[i]->GetYaxis()->SetRangeUser(-20,20);
		kicks_fit_exp17_4_20_res[i]->GetYaxis()->SetRangeUser(-20,20);
	}


	TCanvas* can_fits1 = new TCanvas("can_fits1","",800,800);
	can_fits1->Divide(2,4);
	TCanvas* can_fits1_res = new TCanvas("can_fits1_res","",800,800);
	can_fits1_res->Divide(2,4);

	TCanvas* can_fits2 = new TCanvas("can_fits2","",800,800);
	can_fits2->Divide(2,4);
	TCanvas* can_fits2_res = new TCanvas("can_fits2_res","",800,800);
	can_fits2_res->Divide(2,4);

	TCanvas* can_fits3 = new TCanvas("can_fits3","",800,800);
	can_fits3->Divide(2,4);
	TCanvas* can_fits3_res = new TCanvas("can_fits3_res","",800,800);
	can_fits3_res->Divide(2,4);

	TCanvas* can_fits4 = new TCanvas("can_fits4","",800,800);
	can_fits4->Divide(2,4);
	TCanvas* can_fits4_res = new TCanvas("can_fits4_res","",800,800);
	can_fits4_res->Divide(2,4);

	TGraphErrors* g_off = new TGraphErrors();
	TGraphErrors* g_amp = new TGraphErrors();
	TGraphErrors* g_tau = new TGraphErrors();
	g_off->SetTitle("Offset;Kick;Offset [mV]");
	g_amp->SetTitle("Amplitude;Kick;Amplitude [mV]");
	g_tau->SetTitle("Lifetime;Kick;Lifetime [#mus]");
	g_off->SetMarkerStyle(20);
	g_amp->SetMarkerStyle(20);
	g_tau->SetMarkerStyle(20);
	g_off->SetMarkerColor(kRed);
	g_amp->SetMarkerColor(kRed);
	g_tau->SetMarkerColor(kRed);

	TF1* f_exp = new TF1("f_exp","[2]-[0]*exp(-x/[1])");
	f_exp->SetParameters(10,0.07,0);
	f_exp->SetParNames("A","#tau","offset");

	TF1* f_exp_17 = new TF1("f_exp_17","[2]-[0]*exp(-x/[1])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])");
	f_exp_17->SetParameters(10,0.07,0,2,0.3,17,0);
	f_exp_17->SetParLimits(3,0.5,10);
	f_exp_17->SetParLimits(4,0.1,1.0);
	f_exp_17->SetParLimits(5,15,19);
	f_exp_17->SetParNames("A","#tau","offset","A_{1}","#tau_{1}","f_{1}","#phi_{1}");

	TF1* f_exp_17_4 = new TF1("f_exp_17_4","[2]-[0]*exp(-x/[1])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])+[7]*sin([8]*6.283185307*x+[9])");
	f_exp_17_4->SetParameters(10,0.07,0,2,0.3,17,0,5,4,0);
	f_exp_17_4->SetParLimits(3,0.5,10);
	f_exp_17_4->SetParLimits(4,0.1,1.0);
	f_exp_17_4->SetParLimits(5,15,19);
	f_exp_17_4->SetParLimits(7,0.1,20);
	f_exp_17_4->SetParLimits(8,3,5);
	f_exp_17_4->SetParNames("A","#tau","offset","A_{1}","#tau_{1}","f_{1}","#phi_{1}","A_{2}","f_{2}","#phi_{2}");

	TF1* f_exp_17_4_20 = new TF1("f_exp_17_4_20","[2]-[0]*exp(-x/[1])+[3]*exp(-x/[4])*sin([5]*6.283185307*x+[6])+[7]*sin([8]*6.283185307*x+[9])+([10]+[11]*x)*sin([12]*6.283185307*x+[13])");
	f_exp_17_4_20->SetParameters(10,0.07,0,2,0.3,17,0,5,4,0);
	f_exp_17_4_20->SetParameter(10,10);
	f_exp_17_4_20->SetParameter(11,1);
	f_exp_17_4_20->SetParameter(12,20);
	f_exp_17_4_20->SetParameter(13,0);
	f_exp_17_4_20->SetParLimits(3,0.5,5);
	f_exp_17_4_20->SetParLimits(4,0.1,1.0);
	f_exp_17_4_20->SetParLimits(5,70,120);
	f_exp_17_4_20->SetParLimits(7,0.1,10);
	f_exp_17_4_20->SetParLimits(8,3,5);
	f_exp_17_4_20->SetParLimits(12,18,22);
	f_exp_17_4_20->SetParNames("A","#tau","offset","A_{1}","#tau_{1}","f_{1}","#phi_{1}","A_{2}","f_{2}","#phi_{2}");
	f_exp_17_4_20->SetParName(10,"A_{3}");
	f_exp_17_4_20->SetParName(11,"B_{3}");
	f_exp_17_4_20->SetParName(12,"f_{3}");
	f_exp_17_4_20->SetParName(13,"#phi_{3}");

	for (int i=0; i<8; i++){
		cout<<"Fitting kick "<<i+1<<"\n";
		
		can_fits1->cd(i+1);
		kicks_fit_exp[i]->Draw();
		kicks_fit_exp[i]->GetYaxis()->SetRangeUser(-20,20);

		f_exp->SetParameters(10,0.07,0);
		f_exp->FixParameter(2,0.0);
		TFitResultPtr fit_kick = kicks_fit_exp[i]->Fit(f_exp,"QS","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick->ParError(1));
		if(i==0)cout<<"Amplitude exp: "<<fit_kick->Parameter(0)<<" +- "<<fit_kick->ParError(0)<<"\n";

		kicks_fit_exp_res[i]->Add(f_exp,-1);
		can_fits1_res->cd(i+1);
		kicks_fit_exp_res[i]->Draw();
		kicks_fit_exp_res[i]->GetYaxis()->SetRangeUser(-20,20);
		
		can_fits2->cd(i+1);
		kicks_fit_exp17[i]->Draw();
		kicks_fit_exp17[i]->GetYaxis()->SetRangeUser(-20,20);

		for (int k=0; k<f_exp->GetNpar(); k++){
			f_exp_17->SetParameter(k,f_exp->GetParameter(k));
		}
		TFitResultPtr fit_kick17 = kicks_fit_exp17[i]->Fit(f_exp_17,"QS","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick17->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick17->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick17->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick17->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick17->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick17->ParError(1));
		if(i==0)cout<<"Amplitude exp17: "<<fit_kick17->Parameter(0)<<" +- "<<fit_kick17->ParError(0)<<"\n";

		kicks_fit_exp17_res[i]->Add(f_exp_17,-1);
		can_fits2_res->cd(i+1);
		kicks_fit_exp17_res[i]->Draw();
		kicks_fit_exp17_res[i]->GetYaxis()->SetRangeUser(-20,20);
		
		can_fits3->cd(i+1);
		kicks_fit_exp17_4[i]->Draw();
		kicks_fit_exp17_4[i]->GetYaxis()->SetRangeUser(-20,20);


		for (int k=0; k<f_exp_17->GetNpar(); k++){
			f_exp_17_4->SetParameter(k,f_exp_17->GetParameter(k));
		}
		TFitResultPtr fit_kick17_4 = kicks_fit_exp17_4[i]->Fit(f_exp_17_4,"QS","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick17_4->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick17_4->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick17_4->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick17_4->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick17_4->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick17_4->ParError(1));
		if(i==0)cout<<"Amplitude exp17_4: "<<fit_kick17_4->Parameter(0)<<" +- "<<fit_kick17_4->ParError(0)<<"\n";

		kicks_fit_exp17_4_res[i]->Add(f_exp_17_4,-1);
		can_fits3_res->cd(i+1);
		kicks_fit_exp17_4_res[i]->Draw();
		kicks_fit_exp17_4_res[i]->GetYaxis()->SetRangeUser(-20,20);
		
		can_fits4->cd(i+1);
		kicks_fit_exp17_4_20[i]->Draw();
		kicks_fit_exp17_4_20[i]->GetYaxis()->SetRangeUser(-20,20);


		for (int k=0; k<f_exp_17_4->GetNpar(); k++){
			f_exp_17_4_20->SetParameter(k,f_exp_17_4->GetParameter(k));
		}
		TFitResultPtr fit_kick17_4_20 = kicks_fit_exp17_4_20[i]->Fit(f_exp_17_4_20,"QS","",0.02,0.7);
		g_off->SetPoint(g_off->GetN(),i+1,fit_kick17_4_20->Parameter(2));
		g_off->SetPointError(g_off->GetN()-1,0,fit_kick17_4_20->ParError(2));
		g_amp->SetPoint(g_amp->GetN(),i+1,fit_kick17_4_20->Parameter(0));
		g_amp->SetPointError(g_amp->GetN()-1,0,fit_kick17_4_20->ParError(0));
		g_tau->SetPoint(g_tau->GetN(),i+1,1000*fit_kick17_4_20->Parameter(1));
		g_tau->SetPointError(g_tau->GetN()-1,0,1000*fit_kick17_4_20->ParError(1));
		if(i==0)cout<<"Amplitude exp17_4_20: "<<fit_kick17_4_20->Parameter(0)<<" +- "<<fit_kick17_4_20->ParError(0)<<"\n";

		kicks_fit_exp17_4_20_res[i]->Add(f_exp_17_4_20,-1);
		can_fits4_res->cd(i+1);
		kicks_fit_exp17_4_20_res[i]->Draw();
		kicks_fit_exp17_4_20_res[i]->GetYaxis()->SetRangeUser(-20,20);
	}


	TCanvas* can_fit_values = new TCanvas("can_fit_values","",1500,500);
	can_fit_values->Divide(3,1);
	can_fit_values->cd(1);
	g_off->GetYaxis()->SetRangeUser(-1,1);
	g_off->Draw("APZ");
	gPad->SetGridy();
	can_fit_values->cd(2);
	g_amp->GetYaxis()->SetRangeUser(0,25);
	g_amp->Draw("APZ");
	gPad->SetGridy();
	can_fit_values->cd(3);
	g_tau->GetYaxis()->SetRangeUser(20,100);
	g_tau->Draw("APZ");
	gPad->SetGridy();


	TCanvas* can_fit1_res_fft = new TCanvas("can_fit1_res_fft","",800,800);
	can_fit1_res_fft->Divide(2,2);
	can_fit1_res_fft->cd(1);
	kicks_fit_exp[0]->Draw();
	can_fit1_res_fft->cd(2);
	TH1D* fft_kicks_fit_exp = doFFT(kicks_fit_exp[0],0.1,0.7);
	fft_kicks_fit_exp->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp->Draw("HIST");
	gPad->SetLogx();
	can_fit1_res_fft->cd(3);
	kicks_fit_exp_res[0]->Draw();
	can_fit1_res_fft->cd(4);
	TH1D* fft_kicks_fit_exp_res = doFFT(kicks_fit_exp_res[0],0.1,0.7);
	fft_kicks_fit_exp_res->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp_res->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can_fit2_res_fft = new TCanvas("can_fit2_res_fft","",800,800);
	can_fit2_res_fft->Divide(2,2);
	can_fit2_res_fft->cd(1);
	kicks_fit_exp17[0]->Draw();
	can_fit2_res_fft->cd(2);
	TH1D* fft_kicks_fit_exp17 = doFFT(kicks_fit_exp17[0],0.1,0.7);
	fft_kicks_fit_exp17->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17->Draw("HIST");
	gPad->SetLogx();
	can_fit2_res_fft->cd(3);
	kicks_fit_exp17_res[0]->Draw();
	can_fit2_res_fft->cd(4);
	TH1D* fft_kicks_fit_exp17_res = doFFT(kicks_fit_exp17_res[0],0.1,0.7);
	fft_kicks_fit_exp17_res->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17_res->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can_fit3_res_fft = new TCanvas("can_fit3_res_fft","",800,800);
	can_fit3_res_fft->Divide(2,2);
	can_fit3_res_fft->cd(1);
	kicks_fit_exp17_4[0]->Draw();
	can_fit3_res_fft->cd(2);
	TH1D* fft_kicks_fit_exp17_4 = doFFT(kicks_fit_exp17_4[0],0.1,0.7);
	fft_kicks_fit_exp17_4->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17_4->Draw("HIST");
	gPad->SetLogx();
	can_fit3_res_fft->cd(3);
	kicks_fit_exp17_4_res[0]->Draw();
	can_fit3_res_fft->cd(4);
	TH1D* fft_kicks_fit_exp17_4_res = doFFT(kicks_fit_exp17_4_res[0],0.1,0.7);
	fft_kicks_fit_exp17_4_res->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17_4_res->Draw("HIST");
	gPad->SetLogx();


	TCanvas* can_fit4_res_fft = new TCanvas("can_fit4_res_fft","",800,800);
	can_fit4_res_fft->Divide(2,2);
	can_fit4_res_fft->cd(1);
	kicks_fit_exp17_4_20[0]->Draw();
	can_fit4_res_fft->cd(2);
	TH1D* fft_kicks_fit_exp17_4_20 = doFFT(kicks_fit_exp17_4_20[0],0.1,0.7);
	fft_kicks_fit_exp17_4_20->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17_4_20->Draw("HIST");
	gPad->SetLogx();
	can_fit4_res_fft->cd(3);
	kicks_fit_exp17_4_20_res[0]->Draw();
	can_fit4_res_fft->cd(4);
	TH1D* fft_kicks_fit_exp17_4_20_res = doFFT(kicks_fit_exp17_4_20_res[0],0.1,0.7);
	fft_kicks_fit_exp17_4_20_res->GetYaxis()->SetRangeUser(0,0.03);
	fft_kicks_fit_exp17_4_20_res->Draw("HIST");
	gPad->SetLogx();



	if (saveOutput){
		fout->Write();
		fout->Close();
	}

}
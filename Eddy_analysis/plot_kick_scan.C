#include "analysis_tools.C"



void plot_kick_scan(){

	TFile* f1 = TFile::Open("output_SD_R1_eddy_oct5_H0_Bfield/output.root");
	TFile* f2 = TFile::Open("output_SD_R1_eddy_oct6_H0_Bfield_k583/output.root");
	TFile* f3 = TFile::Open("output_SD_R1_eddy_oct7_H0_nofilter_Bfield_k389/output.root");
	TFile* f4 = TFile::Open("output_SD_R1_eddy_oct7_H0_nofilter_Bfield_k194/output.root");

	TProfile** kicks1 = new TProfile*[8];
	TProfile** kicks2 = new TProfile*[8];
	TProfile** kicks3 = new TProfile*[8];
	TProfile** kicks4 = new TProfile*[8];
	for (int i=0; i<8; i++){
		kicks1[i] = (TProfile*)f1->Get(Form("trace_kick%d",i+1));
		kicks2[i] = (TProfile*)f2->Get(Form("trace_kick%d",i+1));
		kicks3[i] = (TProfile*)f3->Get(Form("trace_kick%d",i+1));
		kicks4[i] = (TProfile*)f4->Get(Form("trace_kick%d",i+1));
	}
	TProfile* kick_average1 = (TProfile*)f1->Get("trace_avg");
	TProfile* kick_average2 = (TProfile*)f2->Get("trace_avg");
	TProfile* kick_average3 = (TProfile*)f3->Get("trace_avg");
	TProfile* kick_average4 = (TProfile*)f4->Get("trace_avg");
	TProfile* trace1 = (TProfile*)f1->Get("trace");
	TProfile* trace2 = (TProfile*)f2->Get("trace");
	TProfile* trace3 = (TProfile*)f3->Get("trace");
	TProfile* trace4 = (TProfile*)f4->Get("trace");

	/*
	for (int bn=1; bn<=trace->GetNbinsX(); bn++){
		if (trace->GetBinContent(bn) < -50){
			trace->SetBinContent(bn,0);
		}
	}
	*/

	gStyle->SetOptStat(0);

	TCanvas* can = new TCanvas("can","",1800,900);
	can->Divide(3,3);
	for (int i=0; i<8; i++){
		can->cd(i+1);
		kicks1[i]->SetLineColor(1);
		kicks2[i]->SetLineColor(2);
		kicks3[i]->SetLineColor(3);
		kicks4[i]->SetLineColor(4);
		kicks1[i]->Draw("HIST L");
		kicks2[i]->Draw("HIST L SAME");
		kicks3[i]->Draw("HIST L SAME");
		kicks4[i]->Draw("HIST L SAME");
		kicks1[i]->GetXaxis()->SetRangeUser(0.0,(i==7?22.0:8.0));
		kicks1[i]->GetYaxis()->SetRangeUser(-20,20);
		gPad->SetGridy();
	}
	can->cd(9);
	kick_average1->SetLineColor(1);
	kick_average2->SetLineColor(2);
	kick_average3->SetLineColor(3);
	kick_average4->SetLineColor(4);
	kick_average1->Draw("HIST L");
	kick_average2->Draw("HIST L SAME");
	kick_average3->Draw("HIST L SAME");
	kick_average4->Draw("HIST L SAME");
	kick_average1->GetXaxis()->SetRangeUser(0,1);
	kick_average1->GetYaxis()->SetRangeUser(-20,20);
	gPad->SetGridy();
	//can->SaveAs("all_kicks.png");

	TCanvas* can2 = new TCanvas("can2","",2500,500);
	trace1->SetMarkerStyle(6);
	trace2->SetMarkerStyle(6);
	trace3->SetMarkerStyle(6);
	trace4->SetMarkerStyle(6);
	trace1->SetLineColor(1);
	trace2->SetLineColor(2);
	trace3->SetLineColor(3);
	trace4->SetLineColor(4);
	trace1->SetMarkerColor(1);
	trace2->SetMarkerColor(2);
	trace3->SetMarkerColor(3);
	trace4->SetMarkerColor(4);
	trace1->GetXaxis()->SetRangeUser(0,100);
	trace1->GetYaxis()->SetRangeUser(-60,60);
	trace1->Draw("P");
	trace2->Draw("P SAME");
	trace3->Draw("P SAME");
	trace4->Draw("P SAME");
	gPad->SetGridy();
	//can2->SaveAs("longtrace.png");


	TCanvas* can_FFT = new TCanvas("can_FFT","",1800,900);
	can_FFT->Divide(3,3);
	for (int i=0; i<9; i++){

		TProfile* this_kick1;
		TProfile* this_kick2;
		TProfile* this_kick3;
		TProfile* this_kick4;
		if (i<8){
			this_kick1 = kicks1[i];
			this_kick2 = kicks2[i];
			this_kick3 = kicks3[i];
			this_kick4 = kicks4[i];
		} else {
			this_kick1 = kick_average1;
			this_kick2 = kick_average2;
			this_kick3 = kick_average3;
			this_kick4 = kick_average4;
		}

		can_FFT->cd(i+1);
		double fft_xmin = 0.1; //ms from kick
		double fft_xmax = (i==7?22.0:8.0); 
    	TH1 *fft_histogram1 = 0;
    	TH1 *fft_histogram2 = 0;
    	TH1 *fft_histogram3 = 0;
    	TH1 *fft_histogram4 = 0;
    	TVirtualFFT::SetTransform(0);
    	
    	TH1F* fftResidualInit1 = SetupFFT(this_kick1, fft_xmin, fft_xmax);
    	fft_histogram1 = fftResidualInit1->FFT(fft_histogram1,"MAG");
    	TH1F* fftResidual1 = RescaleAxis(fft_histogram1, 1./(fft_xmax - fft_xmin));
    	TH1F* fftResidualInit2 = SetupFFT(this_kick2, fft_xmin, fft_xmax);
    	fft_histogram2 = fftResidualInit2->FFT(fft_histogram2,"MAG");
    	TH1F* fftResidual2 = RescaleAxis(fft_histogram2, 1./(fft_xmax - fft_xmin));
    	TH1F* fftResidualInit3 = SetupFFT(this_kick3, fft_xmin, fft_xmax);
    	fft_histogram3 = fftResidualInit3->FFT(fft_histogram3,"MAG");
    	TH1F* fftResidual3 = RescaleAxis(fft_histogram3, 1./(fft_xmax - fft_xmin));
    	TH1F* fftResidualInit4 = SetupFFT(this_kick4, fft_xmin, fft_xmax);
    	fft_histogram4 = fftResidualInit4->FFT(fft_histogram4,"MAG");
    	TH1F* fftResidual4 = RescaleAxis(fft_histogram4, 1./(fft_xmax - fft_xmin));
    	

    	fftResidual1->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	fftResidual1->SetStats(0);
    	fftResidual1->SetName(Form("residualFFT1_%d",i));
    	fftResidual1->Scale(1.0 / fftResidual1->Integral());
    	fftResidual1->GetXaxis()->SetRangeUser(0, fftResidual1->GetXaxis()->GetXmax()/2.);

    	fftResidual2->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	fftResidual2->SetStats(0);
    	fftResidual2->SetName(Form("residualFFT1_%d",i));
    	fftResidual2->Scale(1.0 / fftResidual2->Integral());
    	fftResidual2->GetXaxis()->SetRangeUser(0, fftResidual2->GetXaxis()->GetXmax()/2.);
    	fftResidual3->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	fftResidual3->SetStats(0);
    	fftResidual3->SetName(Form("residualFFT1_%d",i));
    	fftResidual3->Scale(1.0 / fftResidual3->Integral());
    	fftResidual3->GetXaxis()->SetRangeUser(0, fftResidual3->GetXaxis()->GetXmax()/2.);
     	fftResidual4->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	fftResidual4->SetStats(0);
    	fftResidual4->SetName(Form("residualFFT1_%d",i));
    	fftResidual4->Scale(1.0 / fftResidual4->Integral());
    	fftResidual4->GetXaxis()->SetRangeUser(0, fftResidual4->GetXaxis()->GetXmax()/2.);
     	    	
		fftResidual1->SetLineColor(1);
		fftResidual2->SetLineColor(2);
		fftResidual3->SetLineColor(3);
		fftResidual4->SetLineColor(4);

		fftResidual1->Draw("HIST");
    	fftResidual1->GetXaxis()->SetRangeUser(0.2,1000);
    	gPad->SetLogx();
    	fftResidual1->GetXaxis()->SetRangeUser((i==7?0.07:0.2),1000);
    	fftResidual1->GetYaxis()->SetRangeUser(0,0.006);
    	fftResidual2->Draw("HIST SAME");
    	fftResidual3->Draw("HIST SAME");
    	fftResidual4->Draw("HIST SAME");


	}
	//can_FFT->SaveAs("FFT.png");



}
#include "../analysis_tools.C"


void plot_output(TString filename="output.root",TString outfolder=""){

	TFile* f = TFile::Open(filename);

	bool savePlots = false;
	if (outfolder != ""){
		savePlots = true;
	}

	TProfile** kicks = new TProfile*[8];
	TProfile** kicks_SNR0 = new TProfile*[8];
	TProfile** kicks_SNR1 = new TProfile*[8];
	TProfile** kicks_SNR2 = new TProfile*[8];
	for (int i=0; i<8; i++){
		kicks[i] = (TProfile*)f->Get(Form("trace_kick%d",i+1));
		kicks_SNR0[i] = (TProfile*)f->Get(Form("trace_SNR0_kick%d",i+1));
		kicks_SNR1[i] = (TProfile*)f->Get(Form("trace_SNR1_kick%d",i+1));
		kicks_SNR2[i] = (TProfile*)f->Get(Form("trace_SNR2_kick%d",i+1));
	}
	TProfile* trace = (TProfile*)f->Get("trace");
	TProfile* trace_SNR0 = (TProfile*)f->Get("trace_SNR0");
	TProfile* trace_SNR1 = (TProfile*)f->Get("trace_SNR1");
	TProfile* trace_SNR2 = (TProfile*)f->Get("trace_SNR2");

	TProfile* kick8long = (TProfile*)f->Get("trace_kick8long");
	TProfile* kick8long_SNR0 = (TProfile*)f->Get("trace_SNR0_kick8long");
	TProfile* kick8long_SNR1 = (TProfile*)f->Get("trace_SNR1_kick8long");
	TProfile* kick8long_SNR2 = (TProfile*)f->Get("trace_SNR2_kick8long");

	TProfile* kick_average = (TProfile*)f->Get("trace_avg");

	TProfile** g_trace_trend = new TProfile*[10];
	for (int i=0; i<10; i++){
		TString hname = Form("trace_trend_%d_%d",i*10,(i+1)*10);
		g_trace_trend[i] = (TProfile*)f->Get(hname);
	}


	bool haveSNRplots = (trace_SNR0 != nullptr);


	TGraph* g_trend_A = (TGraph*)f->Get("g_trend_A");
	TGraph* g_trend_B = (TGraph*)f->Get("g_trend_B");
	TGraph* g_trend_ABdiff = (TGraph*)f->Get("g_trend_ABdiff");
	TGraphErrors* g_trend_baseline = (TGraphErrors*)f->Get("g_trend_baseline");
	TGraph* g_trend_stddev = (TGraph*)f->Get("g_trend_stddev");
	TGraphErrors* g_trend_blumlein = (TGraphErrors*)f->Get("g_trend_blumlein");
	TGraph* g_trend_SNR = (TGraph*)f->Get("g_trend_SNR");


	TGraph* g_correlation_AB_blumlein = (TGraph*)f->Get("g_correlation_AB_blumlein");
	TGraph* g_correlation_ABdiff_blumlein = (TGraph*)f->Get("g_correlation_ABdiff_blumlein");
	TGraph* g_correlation_AB_SNR = (TGraph*)f->Get("g_correlation_AB_SNR");
	TGraph* g_correlation_ABdiff_SNR = (TGraph*)f->Get("g_correlation_ABdiff_SNR");


	for (int bn=1; bn<=trace->GetNbinsX(); bn++){
		if (trace->GetBinContent(bn) < -50){
			trace->SetBinContent(bn,0);
		}
	}
	for (int i=0; i<10; i++){
		for (int bn=1; bn<=g_trace_trend[i]->GetNbinsX(); bn++){
			if (g_trace_trend[i]->GetBinContent(bn) < -50){
				g_trace_trend[i]->SetBinContent(bn,0);
			}
		}

		for (int bn=1; bn<=g_trace_trend[i]->GetNbinsX(); bn++){
			double val = g_trace_trend[i]->GetBinContent(bn);
			g_trace_trend[i]->SetBinContent(bn,val-i*100);
			g_trace_trend[i]->SetBinEntries(bn,1);
		}
	}

	gStyle->SetOptStat(0);

	TCanvas* can = new TCanvas("can","",1800,900);
	can->Divide(3,3);
	for (int i=0; i<8; i++){
		can->cd(i+1);
		kicks[i]->Draw("HIST L");
		kicks[i]->GetXaxis()->SetRangeUser(0,1);//0.0,(i==7?22.0:8.0));
		kicks[i]->GetYaxis()->SetRangeUser(-20,20);
		gPad->SetGridy();
	}
	can->cd(9);
	kick_average->Draw("HIST L");
	kick_average->GetXaxis()->SetRangeUser(0,1);
	kick_average->GetYaxis()->SetRangeUser(-20,20);
	gPad->SetGridy();
	if(savePlots)can->SaveAs(Form("%s/all_kicks.png",outfolder.Data()));

	TCanvas* can2 = new TCanvas("can2","",2500,500);
	trace->SetMarkerStyle(6);
	trace->Draw("P");
	trace->GetXaxis()->SetRangeUser(0,100);
	trace->GetYaxis()->SetRangeUser(-60,60);
	gPad->SetGridy();
	if(savePlots)can2->SaveAs(Form("%s/full_trace.png",outfolder.Data()));

	TCanvas* can3 = new TCanvas("can3","",2500,500);
	kick_average->SetMarkerStyle(6);
	kick_average->Draw("P");
	kick_average->GetYaxis()->SetRangeUser(-40,20);
	gPad->SetGridy();
	if(savePlots)can3->SaveAs(Form("%s/kick_avg.png",outfolder.Data()));

	gStyle->SetPalette(kRainBow);

	new TCanvas();
	for (int i=0; i<10; i++){
		g_trace_trend[i]->GetYaxis()->SetRangeUser(-1000,100);
		g_trace_trend[i]->Draw("SAME HIST PLC");
	}
	gPad->BuildLegend(0.905,0.14,0.995,0.86);
	g_trace_trend[0]->SetTitle("Trace trend");
	if(savePlots)gPad->SaveAs(Form("%s/trend_traces.png",outfolder.Data()));

	TLine* line = new TLine();

	new TCanvas();
	g_trend_A->SetLineWidth(2);
	g_trend_B->SetLineWidth(2);
	g_trend_ABdiff->SetLineWidth(2);
	g_trend_A->SetLineColor(kBlue);
	g_trend_B->SetLineColor(kRed);
	g_trend_ABdiff->SetLineColor(kBlack);
	g_trend_A->GetYaxis()->SetRangeUser(-1,7);
	g_trend_B->GetYaxis()->SetRangeUser(-1,7);
	g_trend_ABdiff->GetYaxis()->SetRangeUser(-1,7);
	g_trend_A->Draw("APLZ");
	g_trend_B->Draw("PLZ");
	g_trend_ABdiff->Draw("PLZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.3,0.4,0.7,0.6);
	g_trend_A->SetTitle("A B channels trend");
	if(savePlots)gPad->SaveAs(Form("%s/trend_AB_diff.png",outfolder.Data()));

	new TCanvas();
	g_trend_baseline->SetLineWidth(2);
	g_trend_blumlein->SetLineWidth(2);
	g_trend_baseline->SetLineColor(kBlue);
	g_trend_blumlein->SetLineColor(kRed);
	g_trend_baseline->GetYaxis()->SetRangeUser(-30,100);
	g_trend_blumlein->GetYaxis()->SetRangeUser(-30,100);
	g_trend_baseline->Draw("APLZ");
	g_trend_blumlein->Draw("PLZ");
	gPad->SetGridy();
	gPad->BuildLegend(0.2,0.75,0.8,0.85);
	g_trend_baseline->SetTitle("Blumlein and baseline trend");
	line->DrawLine(g_trend_blumlein->GetXaxis()->GetXmin(),0,g_trend_blumlein->GetXaxis()->GetXmax(),0);
	if(savePlots)gPad->SaveAs(Form("%s/trend_blumlein.png",outfolder.Data()));


	new TCanvas();
	g_trend_stddev->SetLineWidth(2);
	g_trend_stddev->Draw("APLZ");
	if(savePlots)gPad->SaveAs(Form("%s/trend_StdDev.png",outfolder.Data()));

	TLine* lineSNR = new TLine();
	new TCanvas();
	g_trend_SNR->SetLineWidth(2);
	g_trend_SNR->Draw("APLZ");
	lineSNR->SetLineColor(kRed);
	lineSNR->SetLineWidth(2);
	lineSNR->SetLineStyle(kDashed);
	lineSNR->DrawLine(g_trend_SNR->GetXaxis()->GetXmin(),10,g_trend_SNR->GetXaxis()->GetXmax(),10);
	lineSNR->DrawLine(g_trend_SNR->GetXaxis()->GetXmin(),15,g_trend_SNR->GetXaxis()->GetXmax(),15);
	if(savePlots)gPad->SaveAs(Form("%s/trend_SNR.png",outfolder.Data()));


	TText* t_text = new TText();
	t_text->SetTextAlign(31);
	TCanvas* can_corr = new TCanvas("can_corr","",900,900);
	can_corr->Divide(2,2);
	can_corr->cd(1);
	g_correlation_AB_blumlein->SetMarkerStyle(20);
	g_correlation_AB_blumlein->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_AB_blumlein->GetCorrelationFactor()*100));
	can_corr->cd(2);
	g_correlation_ABdiff_blumlein->SetMarkerStyle(20);
	g_correlation_ABdiff_blumlein->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_ABdiff_blumlein->GetCorrelationFactor()*100));
	can_corr->cd(3);
	g_correlation_AB_SNR->SetMarkerStyle(20);
	g_correlation_AB_SNR->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_AB_SNR->GetCorrelationFactor()*100));
	can_corr->cd(4);
	g_correlation_ABdiff_SNR->SetMarkerStyle(20);
	g_correlation_ABdiff_SNR->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_ABdiff_SNR->GetCorrelationFactor()*100));
	if(savePlots)can_corr->SaveAs(Form("%s/trend_correlations.png",outfolder.Data()));




	TH1D** h1_fft_kicks = new TH1D* [9];


	TCanvas* can_FFT = new TCanvas("can_FFT","",1800,900);
	can_FFT->Divide(3,3);
	for (int i=0; i<9; i++){

		TProfile* this_kick;
		if (i<8){
			this_kick = kicks[i];
		} else {
			this_kick = kick_average;
		}

		can_FFT->cd(i+1);
		double fft_xmin = 0.1; //ms from kick
		double fft_xmax = 8.0; 
    	TH1 *fft_histogram = 0;
    	TVirtualFFT::SetTransform(0);
    	TH1D* fftResidualInit = SetupFFT(this_kick, fft_xmin, fft_xmax);
    	fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
    	h1_fft_kicks[i] = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
    	h1_fft_kicks[i]->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	h1_fft_kicks[i]->SetStats(0);
    	h1_fft_kicks[i]->SetName(Form("residualFFT_%d",i));
    	h1_fft_kicks[i]->Scale(1.0 / h1_fft_kicks[i]->Integral());
    	h1_fft_kicks[i]->GetXaxis()->SetRangeUser(0, h1_fft_kicks[i]->GetXaxis()->GetXmax()/2.);
    	h1_fft_kicks[i]->Draw("HIST");
    	h1_fft_kicks[i]->GetXaxis()->SetRangeUser(0.2,1000);
    	gPad->SetLogx();
    	h1_fft_kicks[i]->GetXaxis()->SetRangeUser(0.2,1000);
    	h1_fft_kicks[i]->GetYaxis()->SetRangeUser(0,0.008);
	}
	if(savePlots)can_FFT->SaveAs(Form("%s/trace_FFTs.png",outfolder.Data()));


	TCanvas* can_FFT_1_8 = new TCanvas("can_FFT","",900,900);
	can_FFT_1_8->Divide(2,2);
	can_FFT_1_8->cd(1);
	kicks[0]->Draw("HIST L");
	kicks[0]->GetXaxis()->SetRangeUser(0,1);
	kicks[0]->GetYaxis()->SetRangeUser(-30,20);
	gPad->SetGridy();
	can_FFT_1_8->cd(2);
	kicks[7]->Draw("HIST L");
	kicks[7]->GetXaxis()->SetRangeUser(0,1);
	kicks[7]->GetYaxis()->SetRangeUser(-30,20);
	gPad->SetGridy();
	can_FFT_1_8->cd(3);
	h1_fft_kicks[0]->Draw("HIST");
    h1_fft_kicks[0]->GetXaxis()->SetRangeUser(0.2,1000);
    gPad->SetLogx();
    h1_fft_kicks[0]->GetXaxis()->SetRangeUser(0.2,1000);
    h1_fft_kicks[0]->GetYaxis()->SetRangeUser(0,0.008);
	can_FFT_1_8->cd(4);
	h1_fft_kicks[7]->Draw("HIST");
    h1_fft_kicks[7]->GetXaxis()->SetRangeUser(0.2,1000);
    gPad->SetLogx();
    h1_fft_kicks[7]->GetXaxis()->SetRangeUser(0.2,1000);
    h1_fft_kicks[7]->GetYaxis()->SetRangeUser(0,0.008);
	if(savePlots)can_FFT_1_8->SaveAs(Form("%s/trace_FFTs_1_8.png",outfolder.Data()));

/*
	new TCanvas();
	double fft_xmin = 0;
	double fft_xmax = 100; 
    TH1 *fft_histogram = 0;
    TVirtualFFT::SetTransform(0);
    TH1D* fftResidualInit = SetupFFT(trace, fft_xmin, fft_xmax);
    fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
    TH1D* fftResidual = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
    fftResidual->SetTitle("FFT trace;Frequency (kHz);Magnitude [Arb Units]");
    fftResidual->SetStats(0);
    fftResidual->SetName("residualFFT_trace");
    fftResidual->Scale(1.0 / fftResidual->Integral());
    fftResidual->GetXaxis()->SetRangeUser(0, fftResidual->GetXaxis()->GetXmax()/2.);
    fftResidual->Draw("HIST");
    fftResidual->GetXaxis()->SetRangeUser(0.01,1000);
    gPad->SetLogx();
    fftResidual->GetXaxis()->SetRangeUser(0.01,1000);
    fftResidual->GetYaxis()->SetRangeUser(0,0.006);
	//gPad->SaveAs("FFT_trace.png");
*/


	if (haveSNRplots){

		//Adjust kicks
		for (int i=0; i<8; i++){
			for (int bn=1; bn<=kicks_SNR0[i]->GetNbinsX(); bn++){
				//if (kicks_SNR0[i]->GetBinContent(bn) < -50) kicks_SNR0[i]->SetBinContent(bn,0);
				//if (kicks_SNR1[i]->GetBinContent(bn) < -50) kicks_SNR1[i]->SetBinContent(bn,0);
				//if (kicks_SNR2[i]->GetBinContent(bn) < -50) kicks_SNR2[i]->SetBinContent(bn,0);
			}
	
			for (int bn=1; bn<=kicks_SNR0[i]->GetNbinsX(); bn++){
				double val = kicks_SNR1[i]->GetBinContent(bn);
				kicks_SNR1[i]->SetBinContent(bn,val-100);
				kicks_SNR1[i]->SetBinEntries(bn,1);
				val = kicks_SNR2[i]->GetBinContent(bn);
				kicks_SNR2[i]->SetBinContent(bn,val-200);
				kicks_SNR2[i]->SetBinEntries(bn,1);
			}
		}

		//adjust longtraces
		for (int bn=1; bn<=trace_SNR0->GetNbinsX(); bn++){
			//if (trace_SNR0->GetBinContent(bn) < -50) trace_SNR0->SetBinContent(bn,0);
			//if (trace_SNR1->GetBinContent(bn) < -50) trace_SNR1->SetBinContent(bn,0);
			//if (trace_SNR2->GetBinContent(bn) < -50) trace_SNR2->SetBinContent(bn,0);
		}
		for (int bn=1; bn<=trace_SNR0->GetNbinsX(); bn++){
			double val = trace_SNR1->GetBinContent(bn);
			trace_SNR1->SetBinContent(bn,val-100);
			trace_SNR1->SetBinEntries(bn,1);
			val = trace_SNR2->GetBinContent(bn);
			trace_SNR2->SetBinContent(bn,val-200);
			trace_SNR2->SetBinEntries(bn,1);
		}

		//adjust kick8long
		for (int bn=1; bn<=kick8long_SNR0->GetNbinsX(); bn++){
			//if (kick8long_SNR0->GetBinContent(bn) < -50) kick8long_SNR0->SetBinContent(bn,0);
			//if (kick8long_SNR1->GetBinContent(bn) < -50) kick8long_SNR1->SetBinContent(bn,0);
			//if (kick8long_SNR2->GetBinContent(bn) < -50) kick8long_SNR2->SetBinContent(bn,0);
		}
		for (int bn=1; bn<=kick8long_SNR0->GetNbinsX(); bn++){
			double val = kick8long_SNR1->GetBinContent(bn);
			kick8long_SNR1->SetBinContent(bn,val-100);
			kick8long_SNR1->SetBinEntries(bn,1);
			val = kick8long_SNR2->GetBinContent(bn);
			kick8long_SNR2->SetBinContent(bn,val-200);
			kick8long_SNR2->SetBinEntries(bn,1);
		}

		TCanvas* can_trace_SNR = new TCanvas("can_trace_SNR","",1600,800);
		trace_SNR0->GetYaxis()->SetRangeUser(-300,100);
		trace_SNR1->GetYaxis()->SetRangeUser(-300,100);
		trace_SNR2->GetYaxis()->SetRangeUser(-300,100);
		trace_SNR0->Draw("SAME HIST PLC");
		trace_SNR1->Draw("SAME HIST PLC");
		trace_SNR2->Draw("SAME HIST PLC");
		gPad->BuildLegend(0.905,0.14,0.995,0.86);
		trace_SNR0->SetTitle("Full trace");
		if(savePlots)gPad->SaveAs(Form("%s/SNR_traces.png",outfolder.Data()));
	

		TCanvas* can_kick_SNR = new TCanvas("can_kick_SNR","",2600,1000);
		can_kick_SNR->Divide(2,1);
		can_kick_SNR->cd(1);
		kicks_SNR0[0]->GetXaxis()->SetRangeUser(-1,1);
		kicks_SNR1[0]->GetXaxis()->SetRangeUser(-1,1);
		kicks_SNR2[0]->GetXaxis()->SetRangeUser(-1,1);
		kicks_SNR0[0]->GetYaxis()->SetRangeUser(-300,100);
		kicks_SNR1[0]->GetYaxis()->SetRangeUser(-300,100);
		kicks_SNR2[0]->GetYaxis()->SetRangeUser(-300,100);
		kicks_SNR0[0]->Draw("SAME HIST PLC");
		kicks_SNR1[0]->Draw("SAME HIST PLC");
		kicks_SNR2[0]->Draw("SAME HIST PLC");
		gPad->BuildLegend(0.6,0.78,0.88,0.88);
		kicks_SNR0[0]->SetTitle("Kick 1");
		can_kick_SNR->cd(2);
		kick8long_SNR0->GetYaxis()->SetRangeUser(-300,100);
		kick8long_SNR1->GetYaxis()->SetRangeUser(-300,100);
		kick8long_SNR2->GetYaxis()->SetRangeUser(-300,100);
		kick8long_SNR0->Draw("SAME HIST PLC");
		kick8long_SNR1->Draw("SAME HIST PLC");
		kick8long_SNR2->Draw("SAME HIST PLC");
		gPad->BuildLegend(0.6,0.78,0.88,0.88);
		kick8long_SNR0->SetTitle("Kick 8");
		if(savePlots)can_kick_SNR->SaveAs(Form("%s/SNR_kicks.png",outfolder.Data()));


		for (int i=0; i<8; i++){
			for (int bn=1; bn<=kicks_SNR0[i]->GetNbinsX(); bn++){
				double val = kicks_SNR2[i]->GetBinContent(bn);
				kicks_SNR2[i]->SetBinContent(bn,val+100);
				kicks_SNR2[i]->SetBinEntries(bn,1);
			}
		}
		for (int bn=1; bn<=trace_SNR0->GetNbinsX(); bn++){
			double val = trace_SNR2->GetBinContent(bn);
			trace_SNR2->SetBinContent(bn,val+100);
			trace_SNR2->SetBinEntries(bn,1);
		}
		for (int bn=1; bn<=kick8long_SNR0->GetNbinsX(); bn++){
			double val = kick8long_SNR2->GetBinContent(bn);
			kick8long_SNR2->SetBinContent(bn,val+100);
			kick8long_SNR2->SetBinEntries(bn,1);
		}


		TCanvas* can_kick_SNR2 = new TCanvas("can_kick_SNR2","",2600,1000);
		can_kick_SNR2->Divide(2,1);
		can_kick_SNR2->cd(1);
		kicks_SNR0[0]->GetXaxis()->SetRangeUser(-1,1);
		kicks_SNR2[0]->GetXaxis()->SetRangeUser(-1,1);
		kicks_SNR0[0]->GetYaxis()->SetRangeUser(-160,60);
		kicks_SNR2[0]->GetYaxis()->SetRangeUser(-160,60);
		kicks_SNR0[0]->Draw("SAME HIST PLC");
		kicks_SNR2[0]->Draw("SAME HIST PLC");
		gPad->SetGridy();
		gPad->BuildLegend(0.6,0.78,0.88,0.88);
		kicks_SNR0[0]->SetTitle("Kick 1");
		can_kick_SNR2->cd(2);
		kick8long_SNR0->GetYaxis()->SetRangeUser(-160,60);
		kick8long_SNR2->GetYaxis()->SetRangeUser(-160,60);
		kick8long_SNR0->Draw("SAME HIST PLC");
		kick8long_SNR2->Draw("SAME HIST PLC");
		gPad->BuildLegend(0.6,0.78,0.88,0.88);
		gPad->SetGridy();
		kick8long_SNR0->SetTitle("Kick 8");
		if(savePlots)can_kick_SNR2->SaveAs(Form("%s/SNR_kicks_2.png",outfolder.Data()));

	}



	//Save csv files
	if(savePlots){
		ofstream fileout_fulltrace;
		fileout_fulltrace.open(Form("%s/fulltrace.csv",outfolder.Data()),ofstream::trunc);
		fileout_fulltrace<<"Time [ms],Value [mV],Error [mV]\n";
		for (int bn=1; bn<=trace->GetNbinsX(); bn++){
			fileout_fulltrace<<trace->GetBinCenter(bn)<<','<<trace->GetBinContent(bn)<<','<<trace->GetBinError(bn)<<'\n';
		}
		fileout_fulltrace.close();

		for (int i=0; i<8; i++){
			ofstream fileout_kick;
			fileout_kick.open(Form("%s/trace_kick%d.csv",outfolder.Data(),i+1),ofstream::trunc);
			fileout_kick<<"Time [ms],Value [mV],Error [mV]\n";
			for (int bn=1; bn<=kicks[i]->GetNbinsX(); bn++){
				fileout_kick<<kicks[i]->GetBinCenter(bn)<<','<<kicks[i]->GetBinContent(bn)<<','<<kicks[i]->GetBinError(bn)<<'\n';
			}
			fileout_kick.close();
		}

		ofstream fileout_kick8long;
		fileout_kick8long.open(Form("%s/trace_kick8long.csv",outfolder.Data()),ofstream::trunc);
		fileout_kick8long<<"Time [ms],Value [mV],Error [mV]\n";
		for (int bn=1; bn<=kick8long->GetNbinsX(); bn++){
			fileout_kick8long<<kick8long->GetBinCenter(bn)<<','<<kick8long->GetBinContent(bn)<<','<<kick8long->GetBinError(bn)<<'\n';
		}
		fileout_kick8long.close();

		if (haveSNRplots){
			ofstream fileout_fulltrace;
			fileout_fulltrace.open(Form("%s/fulltrace_SNR15.csv",outfolder.Data()),ofstream::trunc);
			fileout_fulltrace<<"Time [ms],Value [mV],Error [mV]\n";
			for (int bn=1; bn<=trace_SNR2->GetNbinsX(); bn++){
				fileout_fulltrace<<trace_SNR2->GetBinCenter(bn)<<','<<trace_SNR2->GetBinContent(bn)<<','<<trace_SNR2->GetBinError(bn)<<'\n';
			}
			fileout_fulltrace.close();
		
			for (int i=0; i<8; i++){
				ofstream fileout_kick;
				fileout_kick.open(Form("%s/trace_SNR15_kick%d.csv",outfolder.Data(),i+1),ofstream::trunc);
				fileout_kick<<"Time [ms],Value [mV],Error [mV]\n";
				for (int bn=1; bn<=kicks_SNR2[i]->GetNbinsX(); bn++){
					fileout_kick<<kicks_SNR2[i]->GetBinCenter(bn)<<','<<kicks_SNR2[i]->GetBinContent(bn)<<','<<kicks_SNR2[i]->GetBinError(bn)<<'\n';
				}
				fileout_kick.close();
			}

			ofstream fileout_kick8long;
			fileout_kick8long.open(Form("%s/trace_SNR15_kick8long.csv",outfolder.Data()),ofstream::trunc);
			fileout_kick8long<<"Time [ms],Value [mV],Error [mV]\n";
			for (int bn=1; bn<=kick8long_SNR2->GetNbinsX(); bn++){
				fileout_kick8long<<kick8long_SNR2->GetBinCenter(bn)<<','<<kick8long_SNR2->GetBinContent(bn)<<','<<kick8long_SNR2->GetBinError(bn)<<'\n';
			}
			fileout_kick8long.close();

		}



	}


}
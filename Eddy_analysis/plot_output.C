
TH1F* SetupFFT(TProfile* h, double xmin, double xmax){
    double timeTot = xmax - xmin;
    double binW = h->GetBinWidth(1);
    int nBins = timeTot / binW;
    TH1F* hout = new TH1F("","",nBins, xmin, xmin + (nBins*binW));
    
    int binCount = 0;
    for (int i(0); i < h->GetXaxis()->GetNbins(); i++){
        if (h->GetBinCenter(i) < xmin ) continue;
        if (binCount > nBins) break;
        binCount++;
        double cont = h->GetBinContent(i);
        double err = h->GetBinError(i);
        hout->SetBinContent(binCount, cont);
        hout->SetBinError(binCount, err);
    }
    return hout;
}

TH1F* RescaleAxis(TH1* input, Double_t Scale) {
    int bins = input->GetNbinsX();
    TAxis* xaxis = input->GetXaxis();
    double* ba = new double[bins+1];
    xaxis->GetLowEdge(ba);
    ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
    for (int i = 0; i < bins+1; i++) {
        ba[i] *= Scale;
    }
    TH1F* out = new TH1F(input->GetName(), input->GetTitle(), bins, ba);
    for (int i = 0; i <= bins; i++) {
        out->SetBinContent(i, input->GetBinContent(i));
        out->SetBinError(i, input->GetBinError(i));
    }
    return out;
}



void plot_output(TString filename="output.root"){

	TFile* f = TFile::Open(filename);

	TProfile** kicks = new TProfile*[8];
	for (int i=0; i<8; i++){
		kicks[i] = (TProfile*)f->Get(Form("trace_kick%d",i+1));
	}
	TProfile* kick_average = (TProfile*)f->Get("trace_avg");
	TProfile* trace = (TProfile*)f->Get("trace");

	TProfile** g_trace_trend = new TProfile*[10];
	for (int i=0; i<10; i++){
		TString hname = Form("trace_trend_%d_%d",i*10,(i+1)*10);
		g_trace_trend[i] = (TProfile*)f->Get(hname);
	}


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
	//can->SaveAs("all_kicks.png");

	TCanvas* can2 = new TCanvas("can2","",2500,500);
	trace->SetMarkerStyle(6);
	trace->Draw("P");
	trace->GetXaxis()->SetRangeUser(0,100);
	trace->GetYaxis()->SetRangeUser(-60,60);
	gPad->SetGridy();
	//can2->SaveAs("longtrace.png");

	TCanvas* can3 = new TCanvas("can3","",2500,500);
	kick_average->SetMarkerStyle(6);
	kick_average->Draw("P");
	kick_average->GetYaxis()->SetRangeUser(-20,20);
	gPad->SetGridy();
	//can3->SaveAs("average.png");

	gStyle->SetPalette(kRainBow);

	new TCanvas();
	for (int i=0; i<10; i++){
		g_trace_trend[i]->GetYaxis()->SetRangeUser(-1000,100);
		g_trace_trend[i]->Draw("SAME HIST PLC");
	}
	gPad->BuildLegend(0.905,0.14,0.995,0.86);



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
	gPad->BuildLegend();

	new TCanvas();
	g_trend_baseline->SetLineWidth(2);
	g_trend_blumlein->SetLineWidth(2);
	g_trend_baseline->SetLineColor(kBlue);
	g_trend_blumlein->SetLineColor(kRed);
	g_trend_baseline->GetYaxis()->SetRangeUser(-30,70);
	g_trend_blumlein->GetYaxis()->SetRangeUser(-30,70);
	g_trend_baseline->Draw("APLZ");
	g_trend_blumlein->Draw("PLZ");
	gPad->BuildLegend();

	new TCanvas();
	g_trend_stddev->SetLineWidth(2);
	g_trend_stddev->Draw("APLZ");
	new TCanvas();
	g_trend_SNR->SetLineWidth(2);
	g_trend_SNR->Draw("APLZ");


	TText* t_text = new TText();
	t_text->SetTextAlign(31);
	TCanvas* g_corr = new TCanvas();
	g_corr->Divide(2,2);
	g_corr->cd(1);
	g_correlation_AB_blumlein->SetMarkerStyle(20);
	g_correlation_AB_blumlein->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_AB_blumlein->GetCorrelationFactor()*100));
	g_corr->cd(2);
	g_correlation_ABdiff_blumlein->SetMarkerStyle(20);
	g_correlation_ABdiff_blumlein->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_ABdiff_blumlein->GetCorrelationFactor()*100));
	g_corr->cd(3);
	g_correlation_AB_SNR->SetMarkerStyle(20);
	g_correlation_AB_SNR->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_AB_SNR->GetCorrelationFactor()*100));
	g_corr->cd(4);
	g_correlation_ABdiff_SNR->SetMarkerStyle(20);
	g_correlation_ABdiff_SNR->Draw("AP");
	t_text->DrawTextNDC(0.89,0.85,Form("Correlation: %.1f%%",g_correlation_ABdiff_SNR->GetCorrelationFactor()*100));






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
    	TH1F* fftResidualInit = SetupFFT(this_kick, fft_xmin, fft_xmax);
    	fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
    	TH1F* fftResidual = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
    	fftResidual->SetTitle(Form("FFT %d;Frequency (kHz);Magnitude [Arb Units]",i+1));
    	fftResidual->SetStats(0);
    	fftResidual->SetName(Form("residualFFT_%d",i));
    	fftResidual->Scale(1.0 / fftResidual->Integral());
    	fftResidual->GetXaxis()->SetRangeUser(0, fftResidual->GetXaxis()->GetXmax()/2.);
    	fftResidual->Draw("HIST");
    	fftResidual->GetXaxis()->SetRangeUser(0.2,1000);
    	gPad->SetLogx();
    	fftResidual->GetXaxis()->SetRangeUser(0.2,1000);
    	fftResidual->GetYaxis()->SetRangeUser(0,0.006);
	}
	//can_FFT->SaveAs("FFT.png");

/*
	new TCanvas();
	double fft_xmin = 0;
	double fft_xmax = 100; 
    TH1 *fft_histogram = 0;
    TVirtualFFT::SetTransform(0);
    TH1F* fftResidualInit = SetupFFT(trace, fft_xmin, fft_xmax);
    fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
    TH1F* fftResidual = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
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

}
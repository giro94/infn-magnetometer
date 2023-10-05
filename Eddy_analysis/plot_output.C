
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



void plot_output(){

	TFile* f = TFile::Open("output.root");

	TProfile** kicks = new TProfile*[8];
	for (int i=0; i<8; i++){
		kicks[i] = (TProfile*)f->Get(Form("trace_kick%d",i+1));
	}
	TProfile* kick_average = (TProfile*)f->Get("trace_avg");
	TProfile* trace = (TProfile*)f->Get("trace");

	gStyle->SetOptStat(0);

	TCanvas* can = new TCanvas("can","",1800,900);
	can->Divide(3,3);
	for (int i=0; i<8; i++){
		can->cd(i+1);
		kicks[i]->Draw("HIST L");
		kicks[i]->GetXaxis()->SetRangeUser(0.0,1.0);
		kicks[i]->GetYaxis()->SetRangeUser(-20,10);
		gPad->SetGridy();
	}
	can->cd(9);
	kick_average->Draw("HIST L");
	kick_average->GetXaxis()->SetRangeUser(0.0,1.0);
	kick_average->GetYaxis()->SetRangeUser(-20,10);
	gPad->SetGridy();
	can->SaveAs("all_kicks.png");

	TCanvas* can2 = new TCanvas("can2","",2500,500);
	trace->SetMarkerStyle(6);
	trace->Draw("P");
	trace->GetXaxis()->SetRangeUser(10,100);
	trace->GetYaxis()->SetRangeUser(-60,60);
	gPad->SetGridy();
	can2->SaveAs("longtrace.png");


	TCanvas* can3 = new TCanvas("can3","",2500,500);
	kick_average->SetMarkerStyle(6);
	kick_average->Draw("P");
	kick_average->GetYaxis()->SetRangeUser(-20,10);
	gPad->SetGridy();
	can3->SaveAs("average.png");

	TFile* fout = new TFile("output_fft.root","recreate");

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
		double fft_xmax = 1.0; 
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
    	fftResidual->GetXaxis()->SetRangeUser(2,1000);
    	gPad->SetLogx();
    	fftResidual->GetXaxis()->SetRangeUser(2,1000);
	}
	can_FFT->SaveAs("FFT.png");


    fout->Write();
    fout->Close();
}
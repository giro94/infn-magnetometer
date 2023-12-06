#include "../analysis_tools.C"

void lowpass(TString input_file){

	TFile* f = TFile::Open(input_file);

	TH1D* trace = ((TProfile*)f->Get("trace_SNR2"))->ProjectionX();
	TH1D** kicks = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks[i] = ((TProfile*)f->Get(Form("trace_SNR2_kick%d",i+1)))->ProjectionX();
	}
	TH1D* kick8long = ((TProfile*)f->Get("trace_SNR2_kick8long"))->ProjectionX();

	for (int bn=1; bn<=trace->GetNbinsX(); bn++){
		if (trace->GetBinContent(bn) < -50){
			trace->SetBinContent(bn,0);
		}
	}

	for (int i=0; i<8; i++){
		for (int bn=1; bn<=kicks[i]->GetNbinsX(); bn++){
			if (kicks[i]->GetBinContent(bn) < -50){
				kicks[i]->SetBinContent(bn,0);
			}
		}
	}

	for (int bn=1; bn<=kick8long->GetNbinsX(); bn++){
		if (kick8long->GetBinContent(bn) < -50){
			kick8long->SetBinContent(bn,0);
		}
	}


	double dt = kicks[0]->GetBinWidth(1);
	double xmin = 0.0;
	double xmax = 2.0;
	int Nbins = (xmax-xmin)/dt;
	TH1D* kick1 = new TH1D("kick1","Kick 1 (original)",Nbins,xmin,xmax);

	for (int bn=1; bn<=kicks[0]->GetNbinsX(); bn++){
		double x = kicks[0]->GetBinCenter(bn);
		if (x >= xmin && x <= xmax){
			kick1->Fill(x,kicks[0]->GetBinContent(bn));
		}
	}


	new TCanvas();
	kick1->Draw("HIST");

	double fft_xmin = xmin;
	double fft_xmax = xmax;

	//Lowpass
	TVirtualFFT::SetTransform(0);

	double timeTot = fft_xmax - fft_xmin;
	int nBins = timeTot / dt;
	TH1F* fftResidualInit = new TH1F("","",nBins, fft_xmin, fft_xmin + (nBins*dt));
	
	int binCount = 0;
	for (int i(0); i < kick1->GetXaxis()->GetNbins(); i++){
		if (kick1->GetBinCenter(i) < fft_xmin ) continue;
		if (binCount > nBins) break;
		binCount++;
		double cont = kick1->GetBinContent(i);
		double err = kick1->GetBinError(i);
		fftResidualInit->SetBinContent(binCount, cont);
		fftResidualInit->SetBinError(binCount, err);
	}

	TH1 *fft_histogram = 0;
	fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();

	TH1 *kick1_re = 0;
	kick1_re = fftResidualInit->FFT(kick1_re, "RE");
	new TCanvas();
	kick1_re->Draw("HIST");
	TH1 *kick1_im = 0;
	kick1_im = fftResidualInit->FFT(kick1_im, "IM");
	new TCanvas();
	kick1_im->Draw("HIST");

	TH1 *kick1_fft = 0;
	kick1_fft = fftResidualInit->FFT(kick1_fft, "MAG");
   	kick1_fft = RescaleAxis(kick1_fft, 1./(fft_xmax - fft_xmin));
	kick1_fft->SetName("kick1_fft");
	kick1_fft->SetTitle("Kick 1 FFT (original)");
	new TCanvas();
	kick1_fft->Draw("HIST");


	double *re_full = new double[Nbins];
	double *im_full = new double[Nbins];
	fft->GetPointsComplex(re_full,im_full);

	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &Nbins, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	TComplex complex_zero = 0;
	double freq_cutoff = 60.0;
	for (int kk=int(freq_cutoff*(fft_xmax - fft_xmin)); kk < Nbins/2 + 1; kk++) {
		fft_back->SetPointComplex(kk, complex_zero);
	}
	fft_back->Transform();


	TH1 *kick1_lowpass = 0;
	kick1_lowpass = TH1::TransformHisto(fft_back,kick1_lowpass,"RE");
   	kick1_lowpass = RescaleAxis(kick1_lowpass, (fft_xmax - fft_xmin)/nBins);
   	kick1_lowpass->Scale(1./nBins);
	kick1_lowpass->SetName("kick1_lowpass");
	kick1_lowpass->SetTitle("Kick 1 (lowpass filtered)");

	new TCanvas();
	kick1_lowpass->Draw("HIST");


	TH1F* fftResidualInit2 = new TH1F("","",nBins, fft_xmin, fft_xmin + (nBins*dt));
	
	binCount = 0;
	for (int i(0); i < kick1_lowpass->GetXaxis()->GetNbins(); i++){
		if (kick1_lowpass->GetBinCenter(i) < fft_xmin ) continue;
		if (binCount > nBins) break;
		binCount++;
		double cont = kick1_lowpass->GetBinContent(i);
		double err = kick1_lowpass->GetBinError(i);
		fftResidualInit2->SetBinContent(binCount, cont);
		fftResidualInit2->SetBinError(binCount, err);
	}
	TH1 *kick1_lowpass_fft = 0;
	kick1_lowpass_fft = fftResidualInit2->FFT(kick1_lowpass_fft, "MAG");
   	kick1_lowpass_fft = RescaleAxis(kick1_lowpass_fft, 1./(fft_xmax - fft_xmin));
	kick1_lowpass_fft->SetName("kick1_lowpass_fft");
	kick1_lowpass_fft->SetTitle("Kick 1 FFT (lowpass filtered)");

	new TCanvas();
	kick1_lowpass_fft->Draw("HIST");


	TCanvas* can = new TCanvas("can","",1200,800);
	can->Divide(2,2);
	can->cd(1);
	kick1->Draw("HIST");
	can->cd(2);
	kick1_fft->Draw("HIST");
	kick1_fft->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
	gPad->SetLogx();
	can->cd(3);
	kick1_lowpass->Draw("HIST");
	can->cd(4);
	kick1_lowpass_fft->Draw("HIST");
	kick1_lowpass_fft->GetXaxis()->SetRangeUser(1./(fft_xmax-fft_xmin),0.5/dt);
	gPad->SetLogx();
	
}
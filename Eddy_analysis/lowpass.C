
TH1D* RescaleAxisLowpass(TH1D* input, double scale, double offset=0) {
    int bins = input->GetNbinsX();
    TAxis* xaxis = input->GetXaxis();
    double* ba = new double[bins+1];
    xaxis->GetLowEdge(ba);
    ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
    for (int i = 0; i < bins+1; i++) {
        ba[i] *= scale;
        ba[i] += offset;
    }
    TH1D* out = new TH1D(Form("%s_rescale",input->GetName()), input->GetTitle(), bins, ba);
    for (int i = 0; i <= bins; i++) {
        out->SetBinContent(i, input->GetBinContent(i));
        out->SetBinError(i, input->GetBinError(i));
    }
    return out;
}

TH1D* lowpass(TH1D* hist, double freq_cutoff, double xmin=-1, double xmax=-1){

	if (xmin == -1) xmin = hist->GetXaxis()->GetXmin();
	if (xmax == -1) xmax = hist->GetXaxis()->GetXmax();

	double binW = hist->GetBinWidth(1);
	double range = xmax-xmin;
	int Nbins = range/binW;

	//Cut input according to args
	TH1D* hist_cut = new TH1D("hist_cut","",Nbins,xmin,xmax);
	for (int bn=1; bn<=hist->GetNbinsX(); bn++){
		double x = hist->GetBinCenter(bn);
		if (x >= xmin && x <= xmax){
			hist_cut->Fill(x,hist->GetBinContent(bn));
		}
	}

	//Perform FFT of input
	TVirtualFFT::SetTransform(0);
	TH1D* hist_fft = 0;
	hist_fft = (TH1D*)hist_cut->FFT(hist_fft,"MAG");

	//Get FFT results in complex space
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	double *re_full = new double[Nbins];
	double *im_full = new double[Nbins];
	fft->GetPointsComplex(re_full,im_full);

	//Perform lowpass filter
	for (int i=int(freq_cutoff*range); i < Nbins/2 + 1; i++) {
		re_full[i] = 0;
		im_full[i] = 0;
	}

	//Perform inverse FFT
	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &Nbins, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();

	//Prepare output
	TH1D *hist_lowpass = 0;
	hist_lowpass = (TH1D*)TH1D::TransformHisto(fft_back,hist_lowpass,"RE");
   	hist_lowpass = RescaleAxisLowpass(hist_lowpass, range/Nbins, xmin);
   	hist_lowpass->Scale(1./Nbins);
	hist_lowpass->SetName(Form("%s_lowpass",hist->GetName()));
	hist_lowpass->SetTitle(Form("%s lowpass",hist->GetTitle()));
   	hist_lowpass->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
   	hist_lowpass->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());

   	//cleanup
   	delete[] re_full;
   	delete[] im_full;
   	delete hist_cut;
   	delete hist_fft;
   	delete fft;
   	delete fft_back;

	return hist_lowpass;	
}
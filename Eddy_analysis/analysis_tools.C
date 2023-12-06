
TH1D* runningAverage(TH1D* input, int width, bool weighted=false){

    TH1D* hist_ra = (TH1D*)input->Clone(Form("%s_runningAverage%d%s",input->GetName(),width,weighted?"W":""));

    int Nbins = hist_ra->GetNbinsX();
    for (int bn=1; bn<=Nbins; bn++){
        int start = bn - width/2;
        int end = bn + width/2;
        if (start < 1) start = 1;
        if (end > Nbins) end = Nbins;
        
        double avg = 0;
        double sumW = 0;
        for (int i=start; i<=end; i++){

            double w = ( weighted ? (1 + width/2 - abs(i-bn)) : 1);
            avg += w*input->GetBinContent(i);
            sumW += w;
        }
        avg /= sumW;

        hist_ra->SetBinContent(bn,avg);
    }
    return hist_ra;
}



TH1D* SetupFFT(TH1D* h_in, double xmin, double xmax){
    double timeTot = xmax - xmin;
    double binW = h_in->GetBinWidth(1);
    int nBins = timeTot / binW;
    TH1D* h_out = new TH1D(Form("%s_fftInit",h_in->GetName()),h_in->GetTitle(),nBins, xmin, xmin + (nBins*binW));
    
    int binCount = 0;
    for (int i=0; i<h_in->GetXaxis()->GetNbins(); i++){
        if (h_in->GetBinCenter(i) < xmin ) continue;
        if (binCount > nBins) break;
        binCount++;
        h_out->SetBinContent(binCount, h_in->GetBinContent(i));
        h_out->SetBinError(binCount, h_in->GetBinError(i));
    }
    return h_out;
}

TH1D* RescaleAxis(TH1* h_in, double scale, double offset=0) {
    int bins = h_in->GetNbinsX();
    TAxis* xaxis = h_in->GetXaxis();
    double* ba = new double[bins+1];
    xaxis->GetLowEdge(ba);
    ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
    for (int i = 0; i < bins+1; i++) {
        ba[i] *= scale;
        ba[i] += offset;
    }
    TH1D* h_out = new TH1D(Form("%s_rescale",h_in->GetName()), h_in->GetTitle(), bins, ba);
    for (int i=0; i<=bins; i++) {
        h_out->SetBinContent(i, h_in->GetBinContent(i));
        h_out->SetBinError(i, h_in->GetBinError(i));
    }
    return h_out;
}

TH1D* doFFT(TH1D* h_in, double xmin, double xmax, TString hname=""){
    double dt = h_in->GetBinWidth(1);
    TH1 *fft_histogram = 0;
    TVirtualFFT::SetTransform(0);
    TH1D* h_in_fftInit = SetupFFT(h_in, xmin, xmax);
    fft_histogram = h_in_fftInit->FFT(fft_histogram,"MAG");
    TH1D* h_out = RescaleAxis(fft_histogram, 1./(xmax - xmin));
    h_out->SetTitle(Form("%s FFT;Frequency (kHz);Magnitude [Arb Units]",h_in->GetTitle()));
    h_out->SetStats(0);
    h_out->SetName(hname=="" ? Form("%s_FFT",h_in->GetName()) : hname);
    h_out->Scale(1.0 / h_out->Integral());
    h_out->GetXaxis()->SetRangeUser(1./(xmax-xmin),0.5/dt);
    delete fft_histogram;
    delete h_in_fftInit;
    return h_out;
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
    hist_lowpass = RescaleAxis(hist_lowpass, range/Nbins, xmin);
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




#include <dirent.h>

vector<TString> getListOfFiles(TString folder){
    vector<TString> files;

    cout<<"Reading folder "<<folder<<"\n";
    DIR* dir = opendir(folder.Data());
    if (!dir) {
        cout<<"Cannot open folder "<<folder<<"!\n";
        throw std::runtime_error{"Cannot open folder"};
    }

    struct dirent* dirfile;
    while((dirfile = readdir(dir)) != NULL){
        if (dirfile->d_type != DT_REG) continue;
        TString fname = dirfile->d_name;
        if (!fname.EndsWith("csv")) continue;
        files.push_back(fname);
    }
    sort(files.begin(),files.end());

    return files;
}

TDatime getFileTime(TString filename){
    int idx = filename.Index(TRegexp("_[0-9]"));
    if (idx >= 0){
        filename.Remove(0,idx+1);
    }
    filename.ReplaceAll(".csv","");
    std::tm tm{};
    std::istringstream iss(filename.Data());
    iss >> std::get_time(&tm, "%m_%d_%Y %H_%M_%S");
    if (iss.fail()) {
        cout<<filename<<"\n";
        cout<<"Can't extract time from file name "<<filename<<"!\n";
        throw std::runtime_error{"failed to parse time string"};
    }
    return TDatime(tm.tm_year+1900,tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
}

pair<int,map<TString,int>> getFileLengthAndHeaders(TString filepath){
    ifstream file;
    file.open(filepath.Data());
    if (!file.is_open()){
        cout<<"cannot open "<<filepath<<"\n";
        throw std::runtime_error{"Cannot open file"};
    }
    char buff[256];
    file.getline(buff,256); //Header, names
    stringstream ss(buff);
    
    //Read headers
    int Nvars = 0;
    map<TString,int> map_varnames;
    while(ss.good()){
        string varname;
        getline(ss,varname,',');
        varname.erase(remove(varname.begin(),varname.end(),'\r'),varname.end());

        for (string s : {" (2)","_Eddy","_rampup","_HWPscan"}){
            size_t idx = varname.find(s);
            if (idx != string::npos){
                varname.erase(idx, s.size());
            }
        }
        cout<<"Found var \""<<varname<<"\"\n";
        map_varnames[varname] = Nvars;
        Nvars++;
    }
    file.getline(buff,256);
    file.getline(buff,256);
    
    //Count lines
    int Nlines = 0;
    while (!file.eof()){
        file.getline(buff,256);
        if (file.eof()) break;
        Nlines++;
    }

    file.close();
    return make_pair(Nlines,map_varnames);
}

vector<vector<double>> readFileTraces(TString filepath, int Nvars){
    ifstream file;
    file.open(filepath.Data());
    if (!file.is_open()){
        cout<<"cannot open "<<filepath<<"\n";
        throw std::runtime_error{"Cannot open file"};
    }
    char buff[256];
    file.getline(buff,256); //Header, names
    file.getline(buff,256); //Header, units
    file.getline(buff,256); //Header, blank row
    char comma;
    vector<vector<double>> vectors;
    vectors.resize(Nvars);
    while (!file.eof()){
        double var=0;
        for (int i=0; i<Nvars; i++){
            if (i>0) file>>comma;
            file>>var;
            if (file.eof()) break;
            vectors[i].push_back(var);
        }
    }
    file.close();
    return vectors;
}

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


TH1D* runningAverage_5_10_15(TH1D* input){
    return runningAverage(runningAverage(runningAverage(input,5,true),10,true),15,true);
}

void cleanTrace(TH1D* trace, double threshold=-50, double setvalue=0){
    for (int bn=1; bn<=trace->GetNbinsX(); bn++){
        if (threshold < 0 && trace->GetBinContent(bn) < threshold){
            trace->SetBinContent(bn,setvalue);
        } else if (threshold > 0 && trace->GetBinContent(bn) > threshold){
            trace->SetBinContent(bn,setvalue);
        } 
    }
    return;
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




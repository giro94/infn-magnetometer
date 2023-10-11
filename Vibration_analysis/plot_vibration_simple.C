#include <dirent.h>

TH1F* SetupFFT(TProfile* h, double xmin, double xmax){
	double timeTot = xmax - xmin;
	double binW = h->GetBinWidth(1);
	int nBins = timeTot / binW;
	TH1F* hout = new TH1F(Form("%s_FFT",h->GetName()),"",nBins, xmin, xmin + (nBins*binW));
	
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


void plot_vibration_simple(TString folder){

	TProfile* p_traceX = new TProfile("p_traceX","",29906,0,100);
	TProfile* p_traceY = new TProfile("p_traceY","",29906,0,100);

	TGraph* g_traceX = new TGraph();
	TGraph* g_traceY = new TGraph();

	DIR* dir = opendir(folder.Data());
	struct dirent* dirfile;
	cout<<"Reading folder "<<folder<<"\n";
	
	vector<TString> files;
	while((dirfile = readdir(dir)) != NULL){
		if (dirfile->d_type != DT_REG) continue;
		TString fname = dirfile->d_name;
		if (!fname.EndsWith("csv")) continue;
		files.push_back(fname);
	}

	sort(files.begin(),files.end());

	int Nfiles = files.size();
	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		vector<pair<double,double>> traceX;
		vector<pair<double,double>> traceY;

		ifstream file;
		TString filename = Form("%s/%s",folder.Data(),fname.Data());
		file.open(filename.Data());
		if (!file.is_open()){
			cout<<"cannot open "<<filename<<"\n";
			return;
		}
		char buff[256];
		file.getline(buff,256); //Header, names

		stringstream ss(buff);
		int Nvar = 0;
		double Apos = -1;
		double Bpos = -1;
		while(ss.good()){
			string substr;
			getline(ss,substr,',');

			if (substr.find("average(A)") != string::npos){
				Apos = Nvar;
			}
			if (substr.find("average(B)") != string::npos){
				Bpos = Nvar;
			}
			Nvar++;
		}

		file.getline(buff,256); //Header, units
		file.getline(buff,256); //Header, blank row

		double time, A, B, C, avgA, avgB;
		char comma;
		while (!file.eof()){
			double var=0;
			for (int i=0; i<Nvar; i++){
				if (i>0) file>>comma;
				file>>var;

				if (i==0){time = var;}
				if (i==Apos){avgA = var;}
				if (i==Bpos){avgB = var;}
			}

			if (file.eof()) break;

			traceX.push_back(make_pair(time,avgA));
			traceY.push_back(make_pair(time,avgB));
		}
		file.close();

		g_traceX->Set(0);
		g_traceY->Set(0);
		for (auto p : traceX){
			time = p.first;
			avgA = p.second;
			g_traceX->AddPoint(time,avgA);
		}
		for (auto p : traceY){
			time = p.first;
			avgB = p.second;
			g_traceY->AddPoint(time,avgB);
		}

		TFitResultPtr resA = g_traceX->Fit("pol0","QNS","",2.0,12.0);
		if (resA>=0){
			double baseline = resA->Parameter(0);
			for (auto p : traceX){
				time = p.first;
				avgA = p.second;
				p_traceX->Fill(time,avgA-baseline);
			}
		}
		TFitResultPtr resB = g_traceY->Fit("pol0","QNS","",2.0,12.0);
		if (resB>=0){
			double baseline = resB->Parameter(0);
			for (auto p : traceY){
				time = p.first;
				avgB = p.second;
				p_traceY->Fill(time,avgB-baseline);
			}
		}
	}

	gStyle->SetOptStat(0);

/*
	TF1* f_linear = new TF1("f_linear","[0]+[1]*x",0,100);

	f_linear->SetParameters(0,1e-08);
	p_traceX->Fit("f_linear","S","",0,100);
	TProfile* p_subtractX = (TProfile*)p_traceX->Clone("p_subtractX");
	p_subtractX->Reset();
	for (int bn=1; bn<p_subtractX->GetNbinsX(); bn++){
		double xcenter = p_subtractX->GetBinCenter(bn);
		double val = p_subtractX->GetBinContent(bn);
		p_traceX->SetBinContent(xcenter,val-f_linear->Eval(xcenter));
	}
	
	p_traceY->Fit("f_linear","S","",0,100);
	TProfile* p_subtractY = (TProfile*)p_traceY->Clone("p_subtractY");
	p_subtractY->Reset();
	for (int bn=1; bn<p_subtractY->GetNbinsY(); bn++){
		double ycenter = p_subtractY->GetBinCenter(bn);
		double val = p_subtractY->GetBinContent(bn);
		p_traceY->SetBinContent(ycenter,val-f_linear->Eval(ycenter));
	}
*/	

	new TCanvas();
	p_traceX->SetLineColor(kBlue);
	p_traceX->SetTitle("X diff");
	p_traceX->GetXaxis()->SetTitle("Time [ms]");
	p_traceX->GetYaxis()->SetTitle("Trace [mV]");
	p_traceX->SetMarkerStyle(8);
	p_traceY->SetLineColor(kRed);
	p_traceY->SetTitle("Y diff");
	p_traceY->GetXaxis()->SetTitle("Time [ms]");
	p_traceY->GetYaxis()->SetTitle("Trace [mV]");
	p_traceY->SetMarkerStyle(8);
	p_traceX->GetYaxis()->SetRangeUser(-10,10);
	p_traceY->GetYaxis()->SetRangeUser(-10,10);
	p_traceX->Draw("HIST");
	p_traceY->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend();




	double fft_xmin = 0.0;
	double fft_xmax = 100.0;
	TH1 *fft_histogramX = 0;
	TH1F* fftResidualInitX = 0;
	TH1F* fftResidualX = 0;
	TH1 *fft_histogramY = 0;
	TH1F* fftResidualInitY = 0;
	TH1F* fftResidualY = 0;
	TVirtualFFT::SetTransform(0);


	fftResidualInitX = SetupFFT(p_traceX, fft_xmin, fft_xmax);
	fft_histogramX = fftResidualInitX->FFT(fft_histogramX,"MAG");
	fftResidualX = RescaleAxis(fft_histogramX, 1./(fft_xmax - fft_xmin));
	fftResidualX->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
	fftResidualX->SetStats(0);
	fftResidualX->SetName(Form("residualFFT X [%.1f,%.1f] ms",fft_xmin, fft_xmax));
	fftResidualX->SetTitle(Form("FFT X [%.1f,%.1f] ms",fft_xmin, fft_xmax));
	fftResidualX->Scale(1.0 / fftResidualX->Integral());
	fftResidualX->GetXaxis()->SetRangeUser(0, fftResidualX->GetXaxis()->GetXmax()/2.);

	fftResidualInitY = SetupFFT(p_traceY, fft_xmin, fft_xmax);
	fft_histogramY = fftResidualInitY->FFT(fft_histogramY,"MAG");
	fftResidualY = RescaleAxis(fft_histogramY, 1./(fft_xmax - fft_xmin));
	fftResidualY->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
	fftResidualY->SetStats(0);
	fftResidualY->SetName(Form("residualFFT Y [%.1f,%.1f] ms",fft_xmin, fft_xmax));
	fftResidualY->SetTitle(Form("FFT Y [%.1f,%.1f] ms",fft_xmin, fft_xmax));
	fftResidualY->Scale(1.0 / fftResidualY->Integral());
	fftResidualY->GetXaxis()->SetRangeUser(0, fftResidualY->GetXaxis()->GetXmax()/2.);


	new TCanvas();
	fftResidualX->SetLineColor(kBlue);
	fftResidualY->SetLineColor(kRed);
	fftResidualY->Draw("HIST");
	fftResidualX->Draw("HIST SAME");
	gPad->SetLogx();
	fftResidualX->GetXaxis()->SetRangeUser(0.015,250);
	fftResidualY->GetXaxis()->SetRangeUser(0.015,250);
	fftResidualX->GetYaxis()->SetRangeUser(0,0.1);
	fftResidualY->GetYaxis()->SetRangeUser(0,0.1);
	gPad->SetGridx();
	gPad->BuildLegend();

}

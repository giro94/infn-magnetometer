#include "../analysis_tools.C"

void plot_vibration_simple(TString folder, TString outfolder=""){


	bool savePlots = false;
	if (outfolder != ""){
		savePlots = true;
	}

	TProfile* p_traceX = new TProfile("p_traceX","",29906,0,100);
	TProfile* p_traceY = new TProfile("p_traceY","",29906,0,100);
	TH1D* p_traceR = new TH1D("p_traceR","",29906,0,100);
	TProfile* p_traceSum = new TProfile("p_traceSum","",29906,0,100);

	TGraph* g_traceX = new TGraph();
	TGraph* g_traceY = new TGraph();

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();


	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		TDatime datetime = getFileTime(fname);
		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> traceX = traces[map_varnames["average(A)"]];
		vector<double> traceY = traces[map_varnames["average(B)"]];
		vector<double> traceSum = traces[map_varnames["Channel C"]];

		g_traceX->Set(0);
		g_traceY->Set(0);
		for (int i=0; i<trace_time.size(); i++){
			g_traceX->SetPoint(i,trace_time[i],traceX[i]);
			g_traceY->SetPoint(i,trace_time[i],traceY[i]);

			p_traceSum->Fill(trace_time[i],traceSum[i]);
		}


		TFitResultPtr resA = g_traceX->Fit("pol0","QNS","",0.0,5.0);
		if (resA>=0){
			double baseline = resA->Parameter(0);
			for (int i=0; i<trace_time.size(); i++){
				p_traceX->Fill(trace_time[i],traceX[i]-baseline);
			}
		}
		TFitResultPtr resB = g_traceY->Fit("pol0","QNS","",0.0,5.0);
		if (resB>=0){
			double baseline = resB->Parameter(0);
			for (int i=0; i<trace_time.size(); i++){
				p_traceY->Fill(trace_time[i],traceY[i]-baseline);
			}
		}
	}

	//adjust for wrong scale
	if (p_traceX->GetMaximum() < 0.1) p_traceX->Scale(1000.);
	if (p_traceY->GetMaximum() < 0.1) p_traceY->Scale(1000.);

	//Calculate total displacement
	for (int bn=1; bn<=p_traceX->GetNbinsX(); bn++){
		double R = sqrt(p_traceX->GetBinContent(bn)*p_traceX->GetBinContent(bn)+p_traceY->GetBinContent(bn)*p_traceY->GetBinContent(bn));
		p_traceR->SetBinContent(bn,p_traceX->GetBinCenter(bn),R);
	}

	//Find kick position
	TSpectrum* spectrum = new TSpectrum();
	spectrum->StaticSearch(p_traceY,2,"nodraw");
	TList *functions = p_traceY->GetListOfFunctions();
	TPolyMarker *peaks = (TPolyMarker*)functions->FindObject("TPolyMarker");
	double* peaks_x = peaks->GetX();
	vector<double> sorted_peaks;
	for (int i=0; i<peaks->GetN(); i++){
		sorted_peaks.push_back(peaks_x[i]);
	}
	sort(sorted_peaks.begin(),sorted_peaks.end());
	double firstkick = sorted_peaks[0];
	double kick_dt = 10;

	if (abs(firstkick-5)>1 && abs(firstkick-10)>1) firstkick = 1e9;

	TLine* line_kicks = new TLine();
	line_kicks->SetLineWidth(2);
	line_kicks->SetLineStyle(kDashed);

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

	double ymin = -5;
	double ymax = 5;

	new TCanvas("","",1800,600);
	p_traceX->SetLineColor(kBlue);
	p_traceX->SetTitle("X displacement");
	p_traceX->GetXaxis()->SetTitle("Time [ms]");
	p_traceX->GetYaxis()->SetTitle("Trace [mV]");
	p_traceX->SetMarkerStyle(8);
	p_traceY->SetLineColor(kRed);
	p_traceY->SetTitle("Y displacement");
	p_traceY->GetXaxis()->SetTitle("Time [ms]");
	p_traceY->GetYaxis()->SetTitle("Trace [mV]");
	p_traceY->SetMarkerStyle(8);
	p_traceX->GetYaxis()->SetRangeUser(ymin,ymax);
	p_traceY->GetYaxis()->SetRangeUser(ymin,ymax);
	p_traceX->Draw("HIST");
	p_traceY->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.75,0.75,0.88,0.88);
	for (int i=0; i<8; i++){
		line_kicks->DrawLine(firstkick+i*kick_dt,ymin,firstkick+i*kick_dt,ymax);
	}
	if(savePlots)gPad->SaveAs(Form("%s/vibrations.png",outfolder.Data()));

	double fft_xmin = 0.0;
	double fft_xmax = 100.0;
	TH1D* fftResidualX = doFFT(p_traceX, fft_xmin, fft_xmax);
	TH1D* fftResidualY = doFFT(p_traceY, fft_xmin, fft_xmax);

	new TCanvas();
	fftResidualX->SetLineColor(kBlue);
	fftResidualY->SetLineColor(kRed);
	fftResidualY->Draw("HIST");
	fftResidualX->Draw("HIST SAME");
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->BuildLegend(0.6,0.7,0.88,0.88);
	if(savePlots)gPad->SaveAs(Form("%s/vibrations_FFT.png",outfolder.Data()));



	new TCanvas("","",900,600);
	TH1D* p_traceX_zoom1 = (TH1D*)p_traceX->Clone("p_traceX_zoom");
	TH1D* p_traceY_zoom1 = (TH1D*)p_traceY->Clone("p_traceY_zoom");
	p_traceX_zoom1->GetXaxis()->SetRangeUser(firstkick,firstkick+1);
	p_traceY_zoom1->GetXaxis()->SetRangeUser(firstkick,firstkick+1);
	p_traceX_zoom1->GetYaxis()->SetRangeUser(-0.2,0.2);
	p_traceY_zoom1->GetYaxis()->SetRangeUser(-0.2,0.2);
	p_traceX_zoom1->Draw("HIST");
	p_traceY_zoom1->Draw("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.6,0.75,0.88,0.88);
	if(savePlots)gPad->SaveAs(Form("%s/firstkick.png",outfolder.Data()));


	TH1D* fftResidualXzoom1 = doFFT(p_traceX_zoom1, 5, 6);
	TH1D* fftResidualYzoom1 = doFFT(p_traceY_zoom1, 5, 6);

	new TCanvas();
	fftResidualXzoom1->SetLineColor(kBlue);
	fftResidualYzoom1->SetLineColor(kRed);
	fftResidualYzoom1->Draw("HIST");
	fftResidualXzoom1->Draw("HIST SAME");
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->BuildLegend(0.6,0.7,0.88,0.88);
	if(savePlots)gPad->SaveAs(Form("%s/firstkick_FFT.png",outfolder.Data()));

}

#include "../analysis_tools.C"

void plot_kicks(TString folder, TString output_file=""){

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();


	TGraph** g_kicks = new TGraph* [Nfiles];


	double baseline_fit_start = 0.0;
	double baseline_fit_end = 1.0;
	double peak_fit_width = 0.05;

	TF1* f_baseline = new TF1("f_baseline","[0]");
	TF1* f_peak = new TF1("f_peak","[0]+[2]*(x-[1])*(x-[1])");

	int Nrepetitions = Nfiles/8;
	TGraphErrors** g_peaks = new TGraphErrors* [Nrepetitions];
	for (int i=0; i<Nrepetitions; i++){
		g_peaks[i] = new TGraphErrors();
	}
	TGraphErrors* g_peaks_time = new TGraphErrors();


	for (int fi=0; fi<Nfiles; fi++){
		
        TString fname = files[fi];
		TDatime datetime = getFileTime(fname);
		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];
		vector<double> trace_A = traces[map_varnames["Channel A"]];
		vector<double> trace_B = traces[map_varnames["Channel B"]];
		vector<double> trace_ABdiff = traces[map_varnames["average(B-A)"]];

		g_kicks[fi] = new TGraph();
		for (int i=0; i<trace_time.size(); i++){
			g_kicks[fi]->SetPoint(i, trace_time[i], trace_ABdiff[i]);
		}      


		int N = g_kicks[fi]->GetN();
		double peak_max = 0;
		double peak_x = 0;
		for (int j=0; j<N; j++){
			if (g_kicks[fi]->GetPointY(j) > peak_max){
				peak_max = g_kicks[fi]->GetPointY(j);
				peak_x = g_kicks[fi]->GetPointX(j);
			}
		}

		cout<<"Peak estimate: x="<<peak_x<<", y="<<peak_max<<"\n";

		f_baseline->SetParameter(0,0);
		TFitResultPtr fit_baseline = g_kicks[fi]->Fit("f_baseline","QS+","",baseline_fit_start,baseline_fit_end);
		f_baseline->Draw("SAME");

		f_peak->SetParameters(0.8,peak_x,-1000.0);
		TFitResultPtr fit_peak = g_kicks[fi]->Fit("f_peak","QS+","",peak_x-peak_fit_width,peak_x+peak_fit_width);
		f_peak->Draw("SAME");

		double baseline = fit_baseline->Parameter(0);
		double peak = fit_peak->Parameter(0);
		double baseline_err = fit_baseline->ParError(0);
		double peak_err = fit_peak->ParError(0);
		double amplitude = abs(peak-baseline);
		double amplitude_err = sqrt(peak_err*peak_err+baseline_err*baseline_err);


		int repetition_idx = (fi-fi%8)/8;
		int kick_idx = fi%8;
		g_peaks[repetition_idx]->SetPoint(kick_idx,kick_idx+1,amplitude);
		g_peaks[repetition_idx]->SetPointError(kick_idx,0,amplitude_err);

		g_peaks_time->SetPoint(fi,datetime.Convert(),amplitude);
		g_peaks_time->SetPointError(fi,0,amplitude_err);

    }

    new TCanvas();
    for (int i=0; i<Nfiles; i++){
    	g_kicks[i]->GetXaxis()->SetTitle("Time [#mus]");
       	g_kicks[i]->GetYaxis()->SetTitle("Trace [V]");
    	g_kicks[i]->GetXaxis()->SetRangeUser(0,5);
    	g_kicks[i]->Draw(i==0?"APL":"PL");
    }

    TCanvas* can_amplitudes = new TCanvas("can_amplitudes","",1800,800);
    can_amplitudes->Divide(2,1);
    can_amplitudes->cd(1);
	g_peaks_time->SetMarkerStyle(20);
	g_peaks_time->SetTitle("Kick amplitude");
	g_peaks_time->GetXaxis()->SetTitle("Time");
	g_peaks_time->GetYaxis()->SetTitle("Kick amplitude [V]");
	g_peaks_time->GetXaxis()->SetTimeFormat("%H:%M");
	g_peaks_time->GetXaxis()->SetTimeOffset(0);
	g_peaks_time->GetXaxis()->SetTimeDisplay(1);
	g_peaks_time->GetYaxis()->SetRangeUser(0.6,0.9);
	g_peaks_time->Draw("APLZ");
	gPad->SetGridy();

	can_amplitudes->cd(2);
	for (int i=0; i<Nrepetitions; i++){
		g_peaks[i]->SetTitle("Kick amplitude");
		g_peaks[i]->GetXaxis()->SetTitle("Kick number");
		g_peaks[i]->GetYaxis()->SetTitle("Kick amplitude [V]");
		g_peaks[i]->GetYaxis()->SetRangeUser(0.6,0.9);
		g_peaks[i]->SetMarkerStyle(20);
		g_peaks[i]->SetMarkerColor(i+1);
		g_peaks[i]->Draw(i==0?"APLZ":"PLZ");
	}
	gPad->SetGridy();







}
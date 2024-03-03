#include "analysis_tools.C"

void average_traces(TString folder, TString foutname=""){

	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	bool saveOutput = false;
	if (foutname != ""){
		saveOutput = true;
	}

	vector<int> colors = {1,2,3,4,9,6};

	//Get headers
	TString firstfile = Form("%s/%s",folder.Data(),files[0].Data());
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(firstfile);
	vector<TString> varunits = getHeaderUnits(firstfile);
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	map<int,TString> map_headers;
    for (const auto& pair : map_varnames) {
        map_headers[pair.second] = pair.first;
    }
    
	vector<vector<double>> traces = readFileTraces(firstfile,Nvars);
	float tstart = traces[map_varnames["Time"]].front();
	float tend = traces[map_varnames["Time"]].back();

	TFile* fout;
	if (saveOutput){
		fout = new TFile(foutname,"recreate");
	}

	TProfile** g_traces = new TProfile* [Nvars-1];
	for (int i=0; i<Nvars-1; i++){
		TString hname = Form("trace_%d",i);
		TString htitle = Form("%s;Time %s;Trace %s",map_headers[i+1].Data(),varunits[0].Data(),varunits[i+1].Data());
		g_traces[i] = new TProfile(hname,htitle,Nlines,tstart,tend);
		g_traces[i]->SetLineColor(colors[i]);
	}

	for (int fi=0; fi<Nfiles; fi++){

		TString fname = files[fi];
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";

		TDatime datetime = getFileTime(fname);
		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[map_varnames["Time"]];

		for (int i=0; i<Nvars-1; i++){
			vector<double> trace = traces[i+1];
			for (int j=0; j<trace_time.size(); j++){
				g_traces[i]->Fill(trace_time[j],trace[j]);
			}
		}

	}

	TH1D** g_traces_ra = new TH1D* [Nvars-1];
	for (int i=0; i<Nvars-1; i++){
		g_traces_ra[i] = runningAverage_5_10_15(g_traces[i]->ProjectionX());
		g_traces_ra[i]->SetLineColor(colors[i]);
	}

	for (int i=0; i<Nvars-1; i++){
		new TCanvas();
		g_traces[i]->Draw("HIST");
	}

	for (int i=0; i<Nvars-1; i++){
		new TCanvas();
		g_traces_ra[i]->Draw("HIST");
	}

	if (saveOutput){
		for (int i=0; i<Nvars-1; i++){
			g_traces[i]->Write();
		}

		fout->Write();
		fout->Close();
	}


}

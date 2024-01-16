#include "analysis_tools.C"

void display_traces(TString fname, TString foutname=""){

	cout<<"Reading file \""<<fname<<"\n";
	ifstream file;
	file.open(fname.Data());
	if (!file.is_open()){
		cout<<"cannot open "<<fname<<"\n";
		return;
	}

	bool saveOutput = false;
	if (foutname != ""){
		saveOutput = true;
	}

	vector<int> colors = {1,2,3,4,9,6};

	//Get headers
	pair<int,map<TString,int>> headers = getFileLengthAndHeaders(fname);
	vector<TString> varunits = getHeaderUnits(fname);
	int Nlines = headers.first;
	map<TString,int> map_varnames = headers.second;
	int Nvars = map_varnames.size();

	map<int,TString> map_headers;
    for (const auto& pair : map_varnames) {
        map_headers[pair.second] = pair.first;
    }
    
	//Get traces
	vector<vector<double>> traces = readFileTraces(fname,Nvars);
	vector<double> trace_time = traces[map_varnames["Time"]];

	file.close();

	TFile* fout;
	if (saveOutput){
		fout = new TFile(foutname,"recreate");
	}

	TGraph** g_traces = new TGraph* [Nvars-1];
	for (int i=0; i<Nvars-1; i++){
		g_traces[i] = new TGraph();

		vector<double> trace = traces[i+1];
		for (int j=0; j<trace_time.size(); j++){
			g_traces[i]->SetPoint(g_traces[i]->GetN(),trace_time[j],trace[j]);
		}
	}

	for (int i=0; i<Nvars-1; i++){
		g_traces[i]->SetTitle(map_headers[i+1]);
		g_traces[i]->GetXaxis()->SetTitle(Form("Time %s",varunits[0].Data()));
		g_traces[i]->GetYaxis()->SetTitle(Form("Trace %s",varunits[i+1].Data()));
		g_traces[i]->SetLineColor(colors[i]);
	}


	for (int i=0; i<Nvars-1; i++){
		new TCanvas();
		g_traces[i]->Draw("AL");
	}

	if (saveOutput){
		for (int i=0; i<Nvars-1; i++){
			g_traces[i]->Write();
		}

		fout->Write();
		fout->Close();
	}


}

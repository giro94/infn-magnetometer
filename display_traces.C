void display_traces(TString fname){

	cout<<"Reading file \""<<fname<<"\n";
	ifstream file;
	file.open(fname.Data());
	if (!file.is_open()){
		cout<<"cannot open "<<fname<<"\n";
		return;
	}

	vector<int> colors = {1,2,3,4,9,6};

	char buff[256];
	file.getline(buff,256); //Header, names
	stringstream ss(buff);
	int Nvar = 0;
	vector<TString> varnames;
	while(ss.good()){
		string substr;
		getline(ss,substr,',');
		substr.erase(remove(substr.begin(),substr.end(),'\r'),substr.end());
		varnames.push_back(substr);
		Nvar++;
	}
	cout<<"Found the following columns:\n";
	for (int i=0; i<Nvar; i++){
		cout<<"\""<<varnames[i]<<"\" ";
	}
	cout<<"\n";

	TGraph** g_traces = new TGraph* [Nvar-1];
	for (int i=0; i<Nvar-1; i++){
		g_traces[i] = new TGraph();
	}

	file.getline(buff,256); //Header, blank row
	while (!file.eof()){
		char comma;
		double var=0;
		vector<double> vars;
		for (int i=0; i<Nvar; i++){
			if (i>0) file>>comma;
			file>>var;
			vars.push_back(var);
		}
		if (file.eof()) break;
		
		for (int i=0; i<Nvar-1; i++){
			g_traces[i]->AddPoint(vars[0],vars[i+1]);
		}
	}
	file.close();

	for (int i=0; i<Nvar-1; i++){
		g_traces[i]->SetTitle(varnames[i+1]);
		g_traces[i]->GetXaxis()->SetTitle("Time [ms]");
		g_traces[i]->GetYaxis()->SetTitle("Trace");
		g_traces[i]->SetLineColor(colors[i]);
		g_traces[i]->GetXaxis()->SetRangeUser(g_traces[i]->GetPointX(0),g_traces[i]->GetPointX(g_traces[i]->GetN()-1));
	}



	if (Nvar-1 == 5){
		g_traces[0]->GetYaxis()->SetRangeUser(0,12);
		g_traces[1]->GetYaxis()->SetRangeUser(0,12);
		g_traces[2]->GetYaxis()->SetRangeUser(-2,2);
		g_traces[3]->GetYaxis()->SetRangeUser(-2,2);
		TCanvas* can = new TCanvas("","",1600,900);
		can->Divide(1,2);
		can->cd(1);
		TVirtualPad* pad_up = gPad;
		pad_up->Divide(2,1);
		pad_up->cd(1);
		g_traces[0]->Draw("AL");
		g_traces[1]->Draw("L");
		gPad->BuildLegend();
		pad_up->cd(2);
		g_traces[2]->Draw("AL");
		g_traces[3]->Draw("L");
		gPad->BuildLegend();
		can->cd(2);
		g_traces[4]->Draw("AL");
	} else {
		for (int i=0; i<Nvar-1; i++){
			new TCanvas("","",1800,800);
			g_traces[i]->Draw("AL");
		}
	}



	new TCanvas("","",1800,800);
	for (int i=0; i<Nvar-1; i++){
		g_traces[i]->Draw(i==0?"AL":"L");
	}


}

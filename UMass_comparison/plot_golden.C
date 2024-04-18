void plot_golden(){


	ifstream fout_R0;
	fout_R0.open("INFN_EC_R0_Bon.csv");

	ifstream fout_R1;
	fout_R1.open("INFN_EC_R1_Bon.csv");

	char comma;

	fout_R0.ignore(256,'\n');
	TGraph* g0 = new TGraph();
	while(!fout_R0.eof()){
		double x,y;
		fout_R0>>x>>comma>>y;
		if (fout_R0.eof()) break;
		g0->AddPoint(x,y);
	}
	fout_R0.close();

	fout_R1.ignore(256,'\n');
	TGraph* g1 = new TGraph();
	while(!fout_R1.eof()){
		double x,y;
		fout_R1>>x>>comma>>y;
		if (fout_R1.eof()) break;
		g1->AddPoint(x,y);
	}
	fout_R1.close();





	ifstream fout_kick_R0;
	fout_kick_R0.open("INFN_kick_R0_Bon.csv");

	ifstream fout_kick_R1;
	fout_kick_R1.open("INFN_kick_R1_Bon.csv");


	fout_kick_R0.ignore(256,'\n');
	TGraph* gk0 = new TGraph();
	while(!fout_kick_R0.eof()){
		double x,y;
		fout_kick_R0>>x>>comma>>y;
		if (fout_kick_R0.eof()) break;
		gk0->AddPoint(x,y);
	}
	fout_kick_R0.close();

	fout_kick_R1.ignore(256,'\n');
	TGraph* gk1 = new TGraph();
	while(!fout_kick_R1.eof()){
		double x,y;
		fout_kick_R1>>x>>comma>>y;
		if (fout_kick_R1.eof()) break;
		gk1->AddPoint(x,y);
	}
	fout_kick_R1.close();

	new TCanvas();
	g0->Draw("AL");
	g1->Draw("L");

	new TCanvas();
	gk0->Draw("AL");
	gk1->Draw("L");

}
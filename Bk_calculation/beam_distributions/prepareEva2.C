{

	ifstream f1("kicker_region_noRF_beamProfile_normalized.csv");
	ifstream f2("kicker_region_xRF_beamProfile_normalized.csv");
	ifstream f3("kicker_region_xyRF5_beamProfile_normalized.csv");
	ifstream f4("kicker_region_xyRF6_beamProfile_normalized.csv");

	TFile* fout = new TFile("beam_dists_eva2.root","recreate");
	TH2D* h2_noRF = new TH2D("noRF","noRF beam distribution",180,-45,45,180,-45,45);
	TH2D* h2_xRF = new TH2D("xRF","xRF beam distribution",180,-45,45,180,-45,45);
	TH2D* h2_xyRF5 = new TH2D("xyRF5","xyRF5 beam distribution",180,-45,45,180,-45,45);
	TH2D* h2_xyRF6 = new TH2D("xyRF6","xyRF6 beam distribution",180,-45,45,180,-45,45);
	h2_noRF->GetXaxis()->SetTitle("x [mm]");
	h2_noRF->GetYaxis()->SetTitle("y [mm]");
	h2_xRF->GetXaxis()->SetTitle("x [mm]");
	h2_xRF->GetYaxis()->SetTitle("y [mm]");
	h2_xyRF5->GetXaxis()->SetTitle("x [mm]");
	h2_xyRF5->GetYaxis()->SetTitle("y [mm]");
	h2_xyRF6->GetXaxis()->SetTitle("x [mm]");
	h2_xyRF6->GetYaxis()->SetTitle("y [mm]");

	f1.ignore(256,'\n');
	while (!f1.eof()){
		double x, y, B;
		char comma;
		f1>>x>>comma>>y>>comma>>B;
		if (f1.eof()) break;
		h2_noRF->Fill(x,y,B);
	}
	f2.ignore(256,'\n');
	while (!f2.eof()){
		double x, y, B;
		char comma;
		f2>>x>>comma>>y>>comma>>B;
		if (f2.eof()) break;
		h2_xRF->Fill(x,y,B);
	}
	f3.ignore(256,'\n');
	while (!f3.eof()){
		double x, y, B;
		char comma;
		f3>>x>>comma>>y>>comma>>B;
		if (f3.eof()) break;
		h2_xyRF5->Fill(x,y,B);
	}
	f4.ignore(256,'\n');
	while (!f4.eof()){
		double x, y, B;
		char comma;
		f4>>x>>comma>>y>>comma>>B;
		if (f4.eof()) break;
		h2_xyRF6->Fill(x,y,B);
	}
	
	fout->Write();
	fout->Close();

}
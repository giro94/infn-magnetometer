{

	TFile* f = TFile::Open("beam_dists.root");

	TH2D* h2_noRF = (TH2D*)f->Get("noRF");
	TH2D* h2_xRF = (TH2D*)f->Get("xRF");
	TH2D* h2_xyRF = (TH2D*)f->Get("xyRF");

	ofstream f_noRF, f_xRF, f_xyRF;
	f_noRF.open("beam_noRF.csv");
	f_xRF.open("beam_xRF.csv");
	f_xyRF.open("beam_xyRF.csv");

	double norm_noRF = h2_noRF->Interpolate(0,0);
	double norm_xRF = h2_xRF->Interpolate(0,0);
	double norm_xyRF = h2_xyRF->Interpolate(0,0);

	f_noRF<<"x [mm],y [mm],Beam dist\n";
	f_xRF<<"x [mm],y [mm],Beam dist\n";
	f_xyRF<<"x [mm],y [mm],Beam dist\n";
	for (double y=-44.75; y<=44.75; y+=0.5){
		for (double x=-44.75; x<=44.75; x+=0.5){
			f_noRF<<x<<","<<y<<","<<h2_noRF->Interpolate(x,y) / norm_noRF<<"\n";
			f_xRF<<x<<","<<y<<","<<h2_xRF->Interpolate(x,y) / norm_xRF<<"\n";
			f_xyRF<<x<<","<<y<<","<<h2_xyRF->Interpolate(x,y) / norm_xyRF<<"\n";
		}
	}

}

{

	ifstream file;
	file.open("umass_2025_03_17_0p1.txt");

	double res = 0.1;
	double xmin = -45 -0.5*res;
	double xmax = +45 +0.5*res;
	int nbins = int(round((xmax-xmin)/res));
	cout<<"resolution: "<<res<<" xmin: "<<xmin<<" xmax: "<<xmax<<" nbins: "<<nbins<<"\n";

	TFile* fout = new TFile("UMass_model_0p1.root","recreate");

	TH2D* h2 = new TH2D("h2","UMass model;x [mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax,nbins,xmin,xmax);
	TH1D* h1_x0 = new TH1D("h1_x0","UMass model [x=0 mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax);
	TH1D* h1_x17p5 = new TH1D("h1_x17p5","UMass model [x=17.5 mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax);
	TH1D* h1_y0 = new TH1D("h1_y0","UMass model [y=0 mm];x [mm];By norm [arb.u.]",nbins,xmin,xmax);

	while(!file.eof()){
		double x, y, by;
		file>>x>>y>>by;
		if (file.eof()) break;

		h2->Fill(x,y,by);
		if (abs(x)<0.1*res){
			h1_x0->Fill(y,by);
		}
		if (abs(x-17.5)<0.1*res){
			h1_x17p5->Fill(y,by);
		}
		if (abs(y)<0.1*res){
			h1_y0->Fill(x,by);
		}
	}

	fout->Write();
	fout->Close();



}
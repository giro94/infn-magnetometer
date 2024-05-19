void plot_beam_distributions(){



	ifstream f_run2, f_run3a, f_run3b;
	f_run2.open("Beam_distribution_run2.csv");
	f_run3a.open("Beam_distribution_run3a.csv");
	f_run3b.open("Beam_distribution_run3b.csv");


	f_run2.ignore(256,'\n');
	f_run3a.ignore(256,'\n');
	f_run3b.ignore(256,'\n');

	TH2D* h2_run2 = new TH2D("h2_run2","Beam distribution Run2;x [mm];y [mm]",180,-45,45,180,-45,45);
	TH2D* h2_run3a = new TH2D("h2_run3a","Beam distribution Run3a;x [mm];y [mm]",180,-45,45,180,-45,45);
	TH2D* h2_run3b = new TH2D("h2_run3b","Beam distribution Run3b;x [mm];y [mm]",180,-45,45,180,-45,45);


	char comma;
	for (int i=0; i<180; i++){
		for (int j=0; j<180; j++){
			double x, y, z;
			f_run2>>x>>comma>>y>>comma>>z;
			h2_run2->Fill(x,y,z);
			f_run3a>>x>>comma>>y>>comma>>z;
			h2_run3a->Fill(x,y,z);
			f_run3b>>x>>comma>>y>>comma>>z;
			h2_run3b->Fill(x,y,z);
		}
	}
	f_run2.close();
	f_run3a.close();
	f_run3b.close();


	TH1D* h1_run2 = h2_run2->ProjectionX("h1_run2");
	TH1D* h1_run3a = h2_run3a->ProjectionX("h1_run3a");
	TH1D* h1_run3b = h2_run3b->ProjectionX("h1_run3b");

	double dR = 17.50;
	double w_crystal = 4;
	double L_crystal = 32;
	TGraph* g_crystal1 = new TGraph();
	TGraph* g_crystal2 = new TGraph();
	g_crystal1->SetPoint(0,-w_crystal/2,-L_crystal/2);
	g_crystal1->SetPoint(1,-w_crystal/2,L_crystal/2);
	g_crystal1->SetPoint(2,w_crystal/2,L_crystal/2);
	g_crystal1->SetPoint(3,w_crystal/2,-L_crystal/2);
	g_crystal1->SetPoint(4,-w_crystal/2,-L_crystal/2);

	g_crystal2->SetPoint(0,dR-w_crystal/2,-L_crystal/2);
	g_crystal2->SetPoint(1,dR-w_crystal/2,L_crystal/2);
	g_crystal2->SetPoint(2,dR+w_crystal/2,L_crystal/2);
	g_crystal2->SetPoint(3,dR+w_crystal/2,-L_crystal/2);
	g_crystal2->SetPoint(4,dR-w_crystal/2,-L_crystal/2);


	TLine* l_R0 = new TLine();
	TLine* l_R1 = new TLine();

	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);
	h2_run2->SetContour(256);
	h2_run3a->SetContour(256);
	h2_run3b->SetContour(256);
	TCanvas* can = new TCanvas("can","",1800,600);
	can->Divide(3,1);
	can->cd(1);
	h2_run2->Draw("colz");
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");
	can->cd(2);
	h2_run3a->Draw("colz");
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");
	can->cd(3);
	h2_run3b->Draw("colz");
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");

	TCanvas* canbis = new TCanvas("canbis","",1800,600);
	canbis->Divide(3,1);
	canbis->cd(1);
	h2_run2->Draw("colz");
	canbis->cd(2);
	h2_run3a->Draw("colz");
	canbis->cd(3);
	h2_run3b->Draw("colz");


	TCanvas* can2 = new TCanvas("can2","",1800,600);
	can2->Divide(3,1);
	can2->cd(1);
	h1_run2->Draw("HIST");
	l_R0->DrawLine(0,0,0,h1_run2->GetMaximum());
	l_R1->DrawLine(dR,0,dR,h1_run2->GetMaximum());
	can2->cd(2);
	h1_run3a->Draw("HIST");
	l_R0->DrawLine(0,0,0,h1_run3a->GetMaximum());
	l_R1->DrawLine(dR,0,dR,h1_run3a->GetMaximum());
	can2->cd(3);
	h1_run3b->Draw("HIST");
	l_R0->DrawLine(0,0,0,h1_run3b->GetMaximum());
	l_R1->DrawLine(dR,0,dR,h1_run3b->GetMaximum());





}
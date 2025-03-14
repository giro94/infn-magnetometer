{


	TFile* f = TFile::Open("BeamDistribution.root");

	TH2D* h2 = (TH2D*)f->Get("h2_beam");
	TH1D* h1x = (TH1D*)f->Get("h_beam_x");
	TH1D* h1y = (TH1D*)f->Get("h_beam_y");

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

	g_crystal1->SetLineColor(kBlack);
	g_crystal2->SetLineColor(kBlack);
	g_crystal1->SetLineWidth(2);
	g_crystal2->SetLineWidth(2);

	TLine* l_R0 = new TLine(0,0,0,1e6);
	TLine* l_R1 = new TLine(dR,0,dR,1e6);
	l_R0->SetLineWidth(2);
	l_R1->SetLineWidth(2);

	gStyle->SetPalette(kViridis);
	gStyle->SetOptStat(0);

	h2->SetContour(128);


	new TCanvas("","",1050,1000);
	h2->Draw("col");
	TLine l;
	l.SetLineColor(kRed);
	l.SetLineStyle(kDashed);
	l.SetLineWidth(1);
	double dx = 90./21.;
	for (int i=0; i<22; i++){
		double x = -45.0 + i*dx;
		l.DrawLine(x,-45,x,45);
	}
	for (int i=0; i<22; i++){
		double y = -45.0 + i*dx;
		l.DrawLine(-45,y,45,y);
	}
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");


	new TCanvas("","",1100,1000);
	h1x->Draw("HIST");
	//h1y->Draw("HIST SAME");
	h1x->SetLineWidth(2);
	h1y->SetLineWidth(2);
	h1x->SetLineColor(kBlue);
	h1y->SetLineColor(kRed);
	h1x->GetYaxis()->SetRangeUser(0,1e6);
	h1y->GetYaxis()->SetRangeUser(0,1e6);
	l_R0->Draw();
	l_R1->Draw();
	TLegend* leg = new TLegend(0.7,0.8,0.9,0.9);
	leg->AddEntry(h1x,"X distribution","L");
	//leg->AddEntry(h1y,"Y distribution","L");
	leg->Draw();
	gPad->SetGridx();

}
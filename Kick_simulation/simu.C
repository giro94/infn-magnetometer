{

	TH1D* h1_x0 = new TH1D("h1_x0","By (x=0);y [mm];By",100,-50,50);
	TH1D* h1_xdR = new TH1D("h1_xdR","By (x=dR);y [mm];By",100,-50,50);
	TH1D* h1_y0 = new TH1D("h1_y0","By (y=0);x [mm];By",100,-50,50);


	TH2D* h2_xy = new TH2D("h2_xy","Field from many wires;x [mm];y [mm];By",100,-50,50,100,-50,50);
	TGraph* g_wires = new TGraph();

	double D = 91.4;
	double R = D/2.;
	double maxTheta = 31. * (3.14159/180.);
	double minr = 1.0;
	double Bmax = 20;
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

	bool first = true;
	for (double x=-R; x<R; x+=1){
		for (double y=-R; y<R; y+=1){
			double B=0;
			for (double a=-maxTheta; a<maxTheta; a+=0.001){

				double w1x = -R*cos(a);
				double w1y = R*sin(a);
				double w2x = R*cos(a);
				double w2y = R*sin(a);
				if (first){
					g_wires->AddPoint(w1x,w1y);
					g_wires->AddPoint(w2x,w2y);
				}

				double r1 = sqrt((x-w1x)*(x-w1x) + (y-w1y)*(y-w1y));
				double r2 = sqrt((x-w2x)*(x-w2x) + (y-w2y)*(y-w2y));
				double By1 = 1/r1;
				double By2 = 1/r2;

				if (r1<minr) By1 = 0;
				if (r2<minr) By2 = 0;

				B += 0.1*By1;
				B += 0.1*By2;

			}
			first = false;

			if (x*x+y*y > R*R) continue;
		
			h2_xy->Fill(x,y,B);

			if (abs(x)<0.5){
				h1_x0->Fill(y,B);
			}
			if (abs(x-dR)<0.5){
				h1_xdR->Fill(y,B);
			}
			if (abs(y)<0.5){
				h1_y0->Fill(x,B);
			}

		}
	}

	cout<<"Using "<<g_wires->GetN()<<" wires\n";

	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);
	h2_xy->SetContour(64);

	new TCanvas("","",800,800);
	g_wires->SetMarkerStyle(20);
	g_wires->SetMarkerColor(kBlack);
	g_crystal1->SetLineWidth(2);
	g_crystal2->SetLineWidth(2);
	g_crystal1->SetMarkerColor(kRed);
	g_crystal2->SetMarkerColor(kRed);
	h2_xy->Draw("COLZ");
	g_wires->Draw("P");
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");

	new TCanvas();
	h1_y0->SetLineWidth(2);
	h1_y0->Draw("HIST");

	new TCanvas();
	h1_x0->GetYaxis()->SetRangeUser(3.5,5.5);
	h1_xdR->GetYaxis()->SetRangeUser(3.5,5.5);
	h1_x0->SetLineWidth(2);
	h1_xdR->SetLineWidth(2);
	h1_x0->SetLineColor(kBlue);
	h1_xdR->SetLineColor(kRed);
	h1_x0->Draw("HIST");
	h1_xdR->Draw("HIST SAME");
	gPad->BuildLegend();


	double B0 = h1_y0->Interpolate(0);
	double BdR = h1_y0->Interpolate(dR);
	cout<<"B at x=0: "<<B0<<"\n";
	cout<<"B at x=dR: "<<BdR<<"\n";
	cout<<"Ratio : "<<BdR/B0<<"\n";


	TFitResultPtr res0 = h1_x0->Fit("pol0","S","",-L_crystal/2,L_crystal/2);
	TFitResultPtr resdR = h1_xdR->Fit("pol0","S","",-L_crystal/2,L_crystal/2);

	double By_R0 = res0->Parameter(0);
	double By_R1 = resdR->Parameter(0);
	cout<<"By in R0 crystal (fit) : "<<By_R0<<"\n";
	cout<<"By in R1 crystal (fit) : "<<By_R1<<"\n";
	cout<<"Ratio : "<<By_R1/By_R0<<"\n";

	double By_R0_int = h1_x0->Integral(h1_x0->FindBin(-L_crystal/2),h1_x0->FindBin(L_crystal/2));
	double By_R1_int = h1_xdR->Integral(h1_xdR->FindBin(-L_crystal/2),h1_xdR->FindBin(L_crystal/2));
	cout<<"By in R0 crystal (integral) : "<<By_R0_int<<"\n";
	cout<<"By in R1 crystal (integral) : "<<By_R1_int<<"\n";
	cout<<"Ratio : "<<By_R1_int/By_R0_int<<"\n";

}
{

	ifstream file;
	file.open("umass_2025_03_17_VC_0p5.txt");

	double xmin = -45;
	double xmax = +45;
	double dx = 0.5;
	int nbins = (xmax-xmin)/dx;

	TH2D* h2 = new TH2D("h2","UMass model;x [mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax,nbins,xmin,xmax);
	TH1D* h1_x0 = new TH1D("h1_x0","UMass model [x=0 mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax);
	TH1D* h1_x17p5 = new TH1D("h1_x17p5","UMass model [x=17.5 mm];y [mm];By norm [arb.u.]",nbins,xmin,xmax);
	TH1D* h1_y0 = new TH1D("h1_y0","UMass model [y=0 mm];x [mm];By norm [arb.u.]",nbins,xmin,xmax);

	while(!file.eof()){
		double x, y, by;
		file>>x>>y>>by;
		if (file.eof()) break;

		h2->Fill(x,y,by);
		if (abs(x)<dx+0.1){
			h1_x0->Fill(y,0.5*by);
		}
		if (abs(x-17.5)<dx+0.1){
			h1_x17p5->Fill(y,0.5*by);
		}
		if (abs(y)<dx+0.1){
			h1_y0->Fill(x,0.5*by);
		}
	}

	gStyle->SetOptStat(0);
	gStyle->SetPalette(kTemperatureMap);
	h2->GetZaxis()->SetRangeUser(-5,5);
	h2->SetContour(256);

	new TCanvas("","",1200,1200);
	h2->Draw("colz");

	h1_x0->SetLineWidth(2);
	h1_x17p5->SetLineWidth(2);
	h1_y0->SetLineWidth(2);

	h1_x0->SetLineColor(kBlue);
	h1_x17p5->SetLineColor(kRed);
	h1_y0->SetLineColor(kBlue);

	h1_x0->GetYaxis()->SetRangeUser(-2,2.5);
	h1_y0->GetYaxis()->SetRangeUser(0,30);
	h2->GetZaxis()->SetRangeUser(-5,5);

	TCanvas* can = new TCanvas("can","",2400,800);
	can->Divide(3,1);
	can->cd(1);
	h2->Draw("colz");
	gPad->SetGridx();
	gPad->SetGridy();
	h2->SetTitle("UMass VC model");
	can->cd(2);
	h1_x0->Draw("HIST");
	h1_x17p5->Draw("HIST SAME");
	gPad->BuildLegend(0.3,0.12,0.7,0.22);
	h1_x0->SetTitle("Vertical slices");
	gPad->SetGridx();
	gPad->SetGridy();
	can->cd(3);
	h1_y0->Draw("HIST");
	gPad->BuildLegend(0.3,0.65,0.7,0.75);
	gPad->SetGridx();
	gPad->SetGridy();
	h1_y0->SetTitle("Radial slice");



	//TFile* fout = new TFile("UMass_model.root","recreate");
	//h2->Write();
	//h1_x0->Write();
	//h1_x17p5->Write();
	//h1_y0->Write();
	//fout->Write();
	//fout->Close();
}
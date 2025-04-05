{

	ifstream file;
	file.open("umass_2025_03_17.txt");

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
		if (abs(x)<0.1){
			h1_x0->Fill(y,by);
		}
		if (abs(x-17.5)<0.1){
			h1_x17p5->Fill(y,by);
		}
		if (abs(y)<0.1){
			h1_y0->Fill(x,by);
		}
	}

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(kTemperatureMap);
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

	TCanvas* can = new TCanvas("can","",2400,800);
	can->Divide(3,1);
	can->cd(1);
	h2->Draw("colz");
	gPad->SetGridx();
	gPad->SetGridy();
	can->cd(2);
	h1_x0->Draw("HIST");
	h1_x17p5->Draw("HIST SAME");
	gPad->BuildLegend();
	gPad->SetGridx();
	gPad->SetGridy();
	can->cd(3);
	h1_y0->Draw("HIST");
	gPad->BuildLegend();
	gPad->SetGridx();
	gPad->SetGridy();



	TFile* fout = new TFile("UMass_model.root","recreate");
	h2->Write();
	h1_x0->Write();
	h1_x17p5->Write();
	h1_y0->Write();
	fout->Write();
	fout->Close();
}
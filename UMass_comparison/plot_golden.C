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

	g0->GetXaxis()->SetTitle("Time [ms]");
	g1->GetXaxis()->SetTitle("Time [ms]");
	g0->GetYaxis()->SetTitle("B field [mG]");
	g1->GetYaxis()->SetTitle("B field [mG]");

	g0->SetLineColor(kBlue);
	g1->SetLineColor(kRed);
	g0->SetLineWidth(2);
	g1->SetLineWidth(2);

	g0->GetXaxis()->SetRangeUser(-0.6,2.0);
	g1->GetXaxis()->SetRangeUser(-0.6,2.0);
	g0->GetYaxis()->SetRangeUser(-100,180);
	g1->GetYaxis()->SetRangeUser(-100,180);


	TGraph* g0_zoom = (TGraph*)g0->Clone("g0_zoom");
	TGraph* g1_zoom = (TGraph*)g1->Clone("g1_zoom");

	TCanvas* can = new TCanvas("can","",1200,900);
	g0->Draw("AL");
	g1->Draw("L");
	gPad->SetGridy();
	can->SetHighLightColor(kBlack);
	TPad* pad_zoom = new TPad("pad_zoom","",0.35,0.5,0.85,0.85);
	pad_zoom->SetBorderMode(1);
	pad_zoom->Draw();
	pad_zoom->cd();
	g0_zoom->GetXaxis()->SetRangeUser(0,0.7);
	g1_zoom->GetXaxis()->SetRangeUser(0,0.7);
	g0_zoom->GetYaxis()->SetRangeUser(-60,20);
	g1_zoom->GetYaxis()->SetRangeUser(-60,20);
	g0_zoom->Draw("AL");
	g1_zoom->Draw("L");
	TLegend* leg = new TLegend(0.4,0.2,0.6,0.4);
	leg->AddEntry(g0_zoom,"R0","L");
	leg->AddEntry(g1_zoom,"R1","L");
	leg->Draw();
	TLine* l30 = new TLine(0.03,-60,0.03,20);
	l30->SetLineWidth(2);
	l30->SetLineStyle(kDashed);
	l30->Draw("SAME");
	TLine* l0 = new TLine(0,0,0.7,0);
	l0->SetLineWidth(2);
	l0->SetLineStyle(kDashed);
	l0->Draw("SAME");

}
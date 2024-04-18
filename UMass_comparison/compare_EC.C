void compare_EC(){

	char comma;

	ifstream f_INFN_EC_R0;
	f_INFN_EC_R0.open("INFN_EC_R0_Bon.csv");
	f_INFN_EC_R0.ignore(256,'\n');
	TGraph* g_INFN_EC_R0 = new TGraph();
	while(!f_INFN_EC_R0.eof()){
		double x,y;
		f_INFN_EC_R0>>x>>comma>>y;
		if (f_INFN_EC_R0.eof()) break;
		g_INFN_EC_R0->AddPoint(x,y);
	}
	f_INFN_EC_R0.close();

	ifstream f_INFN_EC_R1;
	f_INFN_EC_R1.open("INFN_EC_R1_Bon.csv");
	f_INFN_EC_R1.ignore(256,'\n');
	TGraph* g_INFN_EC_R1 = new TGraph();
	while(!f_INFN_EC_R1.eof()){
		double x,y;
		f_INFN_EC_R1>>x>>comma>>y;
		if (f_INFN_EC_R1.eof()) break;
		g_INFN_EC_R1->AddPoint(x,y);
	}
	f_INFN_EC_R1.close();

	ifstream f_INFN_EC_R0_raw;
	f_INFN_EC_R0_raw.open("INFN_EC_R0_Bon_raw.csv");
	f_INFN_EC_R0_raw.ignore(256,'\n');
	TGraph* g_INFN_EC_R0_raw = new TGraph();
	while(!f_INFN_EC_R0_raw.eof()){
		double x,y;
		f_INFN_EC_R0_raw>>x>>comma>>y;
		if (f_INFN_EC_R0_raw.eof()) break;
		g_INFN_EC_R0_raw->AddPoint(x,y);
	}
	f_INFN_EC_R0_raw.close();

	ifstream f_UMass_EC_R0;
	f_UMass_EC_R0.open("UMass_EC_R0_Bon.csv");
	f_UMass_EC_R0.ignore(256,'\n');
	TGraph* g_UMass_EC_R0 = new TGraph();
	while(!f_UMass_EC_R0.eof()){
		double x,y;
		f_UMass_EC_R0>>x>>comma>>y;
		if (f_UMass_EC_R0.eof()) break;
		g_UMass_EC_R0->AddPoint(x,y);
	}
	f_UMass_EC_R0.close();

	g_INFN_EC_R0->GetXaxis()->SetTitle("Time [ms]");
	g_UMass_EC_R0->GetXaxis()->SetTitle("Time [ms]");
	g_INFN_EC_R0->GetYaxis()->SetTitle("B field [mG]");
	g_UMass_EC_R0->GetYaxis()->SetTitle("B field [mG]");

	g_INFN_EC_R0->SetLineColor(kBlue);
	g_UMass_EC_R0->SetLineColor(kRed);
	g_INFN_EC_R0->SetLineWidth(2);
	g_UMass_EC_R0->SetLineWidth(2);

	g_INFN_EC_R0->GetXaxis()->SetRangeUser(-0.6,2.0);
	g_UMass_EC_R0->GetXaxis()->SetRangeUser(-0.6,2.0);
	g_INFN_EC_R0->GetYaxis()->SetRangeUser(-60,140);
	g_UMass_EC_R0->GetYaxis()->SetRangeUser(-60,140);

	TGraph* g_INFN_EC_R0_zoom = (TGraph*)g_INFN_EC_R0->Clone("g_INFN_EC_R0_zoom");
	TGraph* g_UMass_EC_R0_zoom = (TGraph*)g_UMass_EC_R0->Clone("g_UMass_EC_R0_zoom");

	TCanvas* can = new TCanvas("can","",1200,900);
	g_INFN_EC_R0->Draw("AL");
	//g_INFN_EC_R1->Draw("L");
	//g_INFN_EC_R0_raw->Draw("L");
	g_UMass_EC_R0->Draw("L");
	gPad->SetGridy();
	can->SetHighLightColor(kBlack);
	TPad* pad_zoom = new TPad("pad_zoom","",0.35,0.5,0.85,0.85);
	pad_zoom->SetBorderMode(1);
	pad_zoom->Draw();
	pad_zoom->cd();
	g_INFN_EC_R0_zoom->GetXaxis()->SetRangeUser(0,0.7);
	g_UMass_EC_R0_zoom->GetXaxis()->SetRangeUser(0,0.7);
	g_INFN_EC_R0_zoom->GetYaxis()->SetRangeUser(-30,20);
	g_UMass_EC_R0_zoom->GetYaxis()->SetRangeUser(-30,20);
	g_INFN_EC_R0_zoom->Draw("AL");
	g_UMass_EC_R0_zoom->Draw("L");
	TLegend* leg = new TLegend(0.4,0.2,0.6,0.4);
	leg->AddEntry(g_INFN_EC_R0_zoom,"INFN","L");
	leg->AddEntry(g_UMass_EC_R0_zoom,"UMass","L");
	leg->Draw();
	TLine* l30 = new TLine(0.03,-30,0.03,20);
	l30->SetLineWidth(2);
	l30->SetLineStyle(kDashed);
	l30->Draw("SAME");

}
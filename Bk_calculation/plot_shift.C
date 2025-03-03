{

	TF1* f0_0 = new TF1("f0_0","[0]",-45,45);
	f0_0->SetParameters(-16.2);

	TF1* f0 = new TF1("f0","[0]+[1]*x*x",-45,45);
	TF1* f1 = new TF1("f1","[0]+[1]*x*x",-45,45);
	TF1* f2 = new TF1("f2","[0]+[1]*x*x",-45,45);


	f0->SetParameters(-16.2, -0.0630);
	f1->SetParameters(-15.873, -0.0817);
	f2->SetParameters(-15.995, -0.0513);


	TGraph* g0 = new TGraph();
	TGraph* g1 = new TGraph();
	TGraph* g2 = new TGraph();
	TGraph* g0_nofit = new TGraph();
	
	g0->SetPoint(0,0,-16.2);
	g0->SetPoint(1,17.5,-35.5);
	g0->SetPoint(2,45,0);
	g1->SetPoint(0,-2,-16.2);
	g1->SetPoint(1,15.5,-35.5);
	g1->SetPoint(2,45,0);
	g2->SetPoint(0,2,-16.2);
	g2->SetPoint(1,19.5,-35.5);
	g2->SetPoint(2,45,0);

	g0_nofit->SetPoint(0,0,-16.2);
	g0_nofit->SetPoint(1,17.5,-35.5);
	g0_nofit->SetPoint(2,45,0);

	g0->SetMarkerStyle(20);
	g0->SetMarkerColor(kBlue);
	g1->SetMarkerStyle(20);
	g1->SetMarkerColor(kRed);
	g2->SetMarkerStyle(20);
	g2->SetMarkerColor(kGreen);

	g0_nofit->SetMarkerStyle(20);
	g0_nofit->SetMarkerColor(kBlue);

	f0_0->SetLineWidth(2);
	f0_0->SetLineColor(kBlue);
	f0->SetLineWidth(2);
	f0->SetLineColor(kBlue);
	f1->SetLineWidth(2);
	f1->SetLineColor(kRed);
	f2->SetLineWidth(2);
	f2->SetLineColor(kGreen);

	g0->GetXaxis()->SetLimits(-45,45);
	g1->GetXaxis()->SetLimits(-45,45);
	g2->GetXaxis()->SetLimits(-45,45);

	new TCanvas();
	f0->Draw("");
	f1->Draw("SAME");
	f2->Draw("SAME");
	g0->Draw("P");
	g1->Draw("P");
	g2->Draw("P");


	TF1* f0_4 = new TF1("f0_4","[0]+[1]*x*x+[2]*x*x*x*x",-45,45);
	TF1* f1_4 = new TF1("f1_4","[0]+[1]*x*x+[2]*x*x*x*x",-45,45);
	TF1* f2_4 = new TF1("f2_4","[0]+[1]*x*x+[2]*x*x*x*x",-45,45);


	new TCanvas();
	g0->Draw("AP");
	g1->Draw("P");
	g2->Draw("P");


	f0_4->SetLineWidth(2);
	f0_4->SetLineColor(kBlue);
	f1_4->SetLineWidth(2);
	f1_4->SetLineColor(kRed);
	f2_4->SetLineWidth(2);
	f2_4->SetLineColor(kGreen);

	f0_4->SetParameters(-16.2,-0.07529,4.11e-5);
	g0->Fit(f0_4);
	f1_4->SetParameters(-16.2,-0.07529,4.11e-5);
	g1->Fit(f1_4);
	f2_4->SetParameters(-16.2,-0.07529,4.11e-5);
	g2->Fit(f2_4);

	new TCanvas();
	g0_nofit->GetXaxis()->SetLimits(-45,45);
	g0_nofit->Draw("AP");
	f0_0->SetLineColor(kBlack);
	f0->SetLineColor(kRed);
	f0_4->SetLineColor(kGreen);
	f0_0->Draw("SAME");
	f0->Draw("SAME");
	f0_4->Draw("SAME");
	f0_0->SetLineColor(kBlack);
	f0->SetLineColor(kRed);
	f0_4->SetLineColor(kGreen);
	g0_nofit->Draw("P");


}
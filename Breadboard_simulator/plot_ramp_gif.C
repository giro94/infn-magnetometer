void plot_ramp_gif(){

	double degtorad = M_PI/180.;

	TImage *img = TImage::Open("magnetometer.png");

	double xfactor = (double)img->GetHeight()/img->GetWidth();

	TArrow* arrow = new TArrow();
	arrow->SetLineWidth(4);
	arrow->SetLineColor(kGreen);
	arrow->SetFillColor(kGreen);
	double arrowsize = 0.01;
	double aL = 0.05;

	double x1 = 0.8;
	double y1 = 0.82;
	double a1 = 0;

	double x2 = 0.68;
	double y2 = 0.82;
	double a2 = M_PI/4;

	double x3 = 0.38;
	double y3 = 0.54;
	double a3 = M_PI/4;

	double xB1 = 0.09;
	double xB2 = 0.15;
	double yBmax1 = 0.60;
	double yBmin1 = 0.29;
	double yBmax2 = 0.61;
	double yBmin2 = 0.30;

	double x4 = 0.12;
	double y4 = 0.30;
	double a4 = M_PI/4;

	double x5 = 0.38;
	double y5 = 0.68;
	double a5 = M_PI/4;

	double x6 = 0.47;
	double y6 = 0.24;
	double a6 = M_PI/4;

	double x7 = 0.56;
	double y7 = 0.24;
	double a7 = M_PI/4;

	double lambda1 = 22.5;//25;
	double lambda2 = 45;
	lambda1 *= degtorad;
	lambda2 *= degtorad;
	double elR = 0.05;
	double el1x = 1.05;
	double el1y = 0.9;
	double el2x = 1.05;
	double el2y = 0.6;
	TEllipse* el1 = new TEllipse(el1x,el1y,xfactor*elR,elR);
	TEllipse* el2 = new TEllipse(el2x,el2y,xfactor*elR,elR);
	TLine* l1 = new TLine(el1x-xfactor*elR*sin(lambda1),el1y-elR*cos(lambda1),el1x+xfactor*elR*sin(lambda1),el1y+elR*cos(lambda1));
	TLine* l2 = new TLine(el2x-xfactor*elR*sin(lambda2),el2y-elR*cos(lambda2),el2x+xfactor*elR*sin(lambda2),el2y+elR*cos(lambda2));
	el1->SetLineWidth(2);
	el2->SetLineWidth(2);
	l1->SetLineWidth(2);
	l2->SetLineWidth(2);
	TLatex* text_l1 = new TLatex();
	TLatex* text_l2 = new TLatex();
	text_l1->SetTextSize(0.03);
	text_l1->SetTextAlign(22);
	text_l2->SetTextSize(0.03);
	text_l2->SetTextAlign(22);

	double B = 0;
	double I = 1;
	double chA = 0;
	double chB = 0;
	double chC = 0;

	TLine* lineA = new TLine(0.60,0.24,0.73,0.18);
	TLine* lineA2 = new TLine(0.73,0.18,0.64,0.40);
	TLine* lineB = new TLine(0.60,0.24,0.79,0.25);
	TLine* lineB2 = new TLine(0.79,0.25,0.67,0.46);
	lineA->SetLineColor(kRed);
	lineA2->SetLineColor(kRed);
	lineB->SetLineColor(kRed);
	lineB2->SetLineColor(kRed);
	double wA = 0;
	double wB = 0;
	double wMax = 20;

	TGraph* g_AB_ramp = new TGraph();
	g_AB_ramp->SetLineWidth(2);
	g_AB_ramp->GetXaxis()->SetTitle("B field [T]");
	g_AB_ramp->GetXaxis()->SetTitleSize(0.08);
	g_AB_ramp->GetYaxis()->SetTitle("Channel difference");
	g_AB_ramp->GetYaxis()->SetTitleSize(0.08);

	TH1D* h1_B = new TH1D("h1_B","",3,0,3);
	h1_B->GetYaxis()->SetTitle("B field [T]");
	h1_B->GetYaxis()->SetTitleSize(0.15);
	h1_B->GetYaxis()->SetLabelSize(0.15);
	h1_B->GetYaxis()->SetLabelOffset(0.04);
	h1_B->GetXaxis()->SetLabelSize(0);
	h1_B->GetXaxis()->SetTickLength(0);
	h1_B->GetYaxis()->SetRangeUser(0,1.5);
	h1_B->SetFillColor(kBlue);
	h1_B->SetLineColor(kBlue);
	
	gStyle->SetOptStat(0);
	TCanvas* can = new TCanvas("can","",img->GetWidth(),img->GetHeight());

	int framecount=0;
	for (B=0; B<=1.45; B+=0.01){
		a1 = 0;
		a2 = 2*lambda1-a1;
		a3 = a2;
		a4 = a3 + B*131*0.032;
		a5 = a3 + B*131*0.032*2;
		a6 = a5;
		a7 = 2*lambda2-a6;
		chA = I*sin(a7)*sin(a7);
		chB = I*cos(a7)*cos(a7);
		chC = chB-chA;
		wA = wMax*chA/I;
		wB = wMax*chB/I;
		g_AB_ramp->AddPoint(B,chC);
		h1_B->SetBinContent(2,B);

		can->Clear();
		img->Draw();
		arrow->DrawArrow(x1-xfactor*aL*sin(a1),y1-aL*cos(a1),x1+xfactor*aL*sin(a1),y1+aL*cos(a1),arrowsize);
		arrow->DrawArrow(x2-xfactor*aL*sin(a2),y2-aL*cos(a2),x2+xfactor*aL*sin(a2),y2+aL*cos(a2),arrowsize);
		arrow->DrawArrow(x3-xfactor*aL*sin(a3),y3-aL*cos(a3),x3+xfactor*aL*sin(a3),y3+aL*cos(a3),arrowsize);
		//arrow->DrawArrow(x4-xfactor*aL*sin(a4),y4-aL*cos(a4),x4+xfactor*aL*sin(a4),y4+aL*cos(a4),arrowsize);
		arrow->DrawArrow(x5-xfactor*aL*sin(a5),y5-aL*cos(a5),x5+xfactor*aL*sin(a5),y5+aL*cos(a5),arrowsize);
		arrow->DrawArrow(x6-xfactor*aL*sin(a6),y6-aL*cos(a6),x6+xfactor*aL*sin(a6),y6+aL*cos(a6),arrowsize);
		arrow->DrawArrow(x7-xfactor*aL*sin(a7),y7-aL*cos(a7),x7+xfactor*aL*sin(a7),y7+aL*cos(a7),arrowsize);
		for (double yB=yBmax1; yB>=yBmin1; yB-=0.1){
			double aB = a3 + ((yB-yBmax1)/(yBmin1-yBmax1))*(a4-a3);
			arrow->DrawArrow(xB1-xfactor*aL*sin(aB),yB-aL*cos(aB),xB1+xfactor*aL*sin(aB),yB+aL*cos(aB),arrowsize);	
		}
		for (double yB=yBmin2; yB<=yBmax2; yB+=0.1){
			double aB = a4 + ((yB-yBmin2)/(yBmax2-yBmin2))*(a5-a4);
			arrow->DrawArrow(xB2-xfactor*aL*sin(aB),yB-aL*cos(aB),xB2+xfactor*aL*sin(aB),yB+aL*cos(aB),arrowsize);	
		}

		lineA->SetLineWidth(wA);
		lineA2->SetLineWidth(wA);
		lineB->SetLineWidth(wB);
		lineB2->SetLineWidth(wB);
		lineA->Draw("SAME");
		lineA2->Draw("SAME");
		lineB->Draw("SAME");
		lineB2->Draw("SAME");
		el1->Draw("SAME");
		el2->Draw("SAME");
		l1->Draw("SAME");
		l2->Draw("SAME");
		text_l1->DrawLatex(el1x,el1y+0.1,Form("HWPin: %.1f#circ",lambda1/degtorad));
		text_l2->DrawLatex(el2x,el2y+0.1,Form("HWPout: %.1f#circ",lambda2/degtorad));
	
		TPad* padB = new TPad("padB","",0,0.3,0.1,0.7);
		padB->Draw();
		padB->cd();
		gPad->SetLeftMargin(0.5);
		h1_B->Draw("HIST F");

		can->cd();
		TPad* pad1 = new TPad("pad1","",0.75,0,1.0,0.4);
		pad1->Draw();
		pad1->cd();
		g_AB_ramp->GetXaxis()->SetLimits(0,1.5);
		g_AB_ramp->GetXaxis()->SetRangeUser(0,1.5);
		g_AB_ramp->GetYaxis()->SetRangeUser(-1.1,1.1);
		gPad->SetBottomMargin(0.2);
		gPad->SetLeftMargin(0.2);
		g_AB_ramp->Draw("APL");
		gPad->SetGridy();

		can->Update();
		can->SaveAs(Form("gif_rampup/%04d.png",framecount++));
		
	}

	framecount=0;
	for (; abs(chC)>0.01; lambda2-=0.15*degtorad){
		a1 = 0;
		a2 = 2*lambda1-a1;
		a3 = a2;
		a4 = a3 + B*131*0.032;
		a5 = a3 + B*131*0.032*2;
		a6 = a5;
		a7 = 2*lambda2-a6;
		chA = I*sin(a7)*sin(a7);
		chB = I*cos(a7)*cos(a7);
		chC = chB-chA;
		wA = wMax*chA/I;
		wB = wMax*chB/I;
		g_AB_ramp->AddPoint(B,chC);
		h1_B->SetBinContent(2,B);

		can->Clear();
		img->Draw();
		arrow->DrawArrow(x1-xfactor*aL*sin(a1),y1-aL*cos(a1),x1+xfactor*aL*sin(a1),y1+aL*cos(a1),arrowsize);
		arrow->DrawArrow(x2-xfactor*aL*sin(a2),y2-aL*cos(a2),x2+xfactor*aL*sin(a2),y2+aL*cos(a2),arrowsize);
		arrow->DrawArrow(x3-xfactor*aL*sin(a3),y3-aL*cos(a3),x3+xfactor*aL*sin(a3),y3+aL*cos(a3),arrowsize);
		arrow->DrawArrow(x4-xfactor*aL*sin(a4),y4-aL*cos(a4),x4+xfactor*aL*sin(a4),y4+aL*cos(a4),arrowsize);
		arrow->DrawArrow(x5-xfactor*aL*sin(a5),y5-aL*cos(a5),x5+xfactor*aL*sin(a5),y5+aL*cos(a5),arrowsize);
		arrow->DrawArrow(x6-xfactor*aL*sin(a6),y6-aL*cos(a6),x6+xfactor*aL*sin(a6),y6+aL*cos(a6),arrowsize);
		arrow->DrawArrow(x7-xfactor*aL*sin(a7),y7-aL*cos(a7),x7+xfactor*aL*sin(a7),y7+aL*cos(a7),arrowsize);

		lineA->SetLineWidth(wA);
		lineA2->SetLineWidth(wA);
		lineB->SetLineWidth(wB);
		lineB2->SetLineWidth(wB);
		lineA->Draw("SAME");
		lineA2->Draw("SAME");
		lineB->Draw("SAME");
		lineB2->Draw("SAME");
		el1->Draw("SAME");
		el2->Draw("SAME");
		l1->Draw("SAME");
		l2->Draw("SAME");
		text_l1->DrawLatex(el1x,el1y+0.1,Form("HWPin: %.1f#circ",lambda1/degtorad));
		text_l2->DrawLatex(el2x,el2y+0.1,Form("HWPout: %.1f#circ",lambda2/degtorad));
	
		TPad* padB = new TPad("padB","",0,0.3,0.1,0.7);
		padB->Draw();
		padB->cd();
		gPad->SetLeftMargin(0.5);
		h1_B->Draw("HIST F");


		can->cd();
		TPad* pad1 = new TPad("pad1","",0.75,0,1.0,0.4);
		pad1->Draw();
		pad1->cd();
		g_AB_ramp->GetXaxis()->SetLimits(0,1.5);
		g_AB_ramp->GetXaxis()->SetRangeUser(0,1.5);
		g_AB_ramp->GetYaxis()->SetRangeUser(-1.1,1.1);
		gPad->SetBottomMargin(0.2);
		gPad->SetLeftMargin(0.2);
		g_AB_ramp->Draw("APL");
		gPad->SetGridy();

		can->Update();
		can->SaveAs(Form("gif_HWPout/%04d.png",framecount++));
		
	}


	ifstream f_kick;
	f_kick.open("../UMass_comparison/INFN_kick_R0_Bon.csv");
	f_kick.ignore(256,'\n');
	vector<double> kick_x;
	vector<double> kick_y;
	while(!f_kick.eof()){
		double x,y;
		char comma;
		f_kick>>x>>comma>>y;
		if (f_kick.eof()) break;
		kick_x.push_back(x);
		kick_y.push_back(y);
	}
	f_kick.close();


	TGraph* g_AB_kick = new TGraph();
	g_AB_kick->SetLineWidth(2);
	g_AB_kick->GetXaxis()->SetTitle("Time [#mus]");
	g_AB_kick->GetXaxis()->SetTitleSize(0.08);
	g_AB_kick->GetYaxis()->SetTitle("Channel difference");
	g_AB_kick->GetYaxis()->SetTitleSize(0.08);

	double B0 = B;
	int round = 1;
	framecount = 0;
	for (int i=0; i<kick_x.size(); i+=5){
		if (kick_x[i]<-1) continue;
		if (kick_x[i]>4) {
			if (round<1){
				i=0;
				round++;
				g_AB_kick->Set(0);
			}
			continue;
		}

		double Bk = 0.03*kick_y[i];
		B = B0 - Bk;
		a1 = 0;
		a2 = 2*lambda1-a1;
		a3 = a2;
		a4 = a3 + B*131*0.032;
		a5 = a3 + B*131*0.032*2;
		a6 = a5;
		a7 = 2*lambda2-a6;
		chA = I*sin(a7)*sin(a7);
		chB = I*cos(a7)*cos(a7);
		chC = chB-chA;
		wA = wMax*chA/I;
		wB = wMax*chB/I;
		g_AB_kick->AddPoint(kick_x[i],kick_y[i]);
		h1_B->SetBinContent(2,B);

		can->Clear();
		img->Draw();
		arrow->DrawArrow(x1-xfactor*aL*sin(a1),y1-aL*cos(a1),x1+xfactor*aL*sin(a1),y1+aL*cos(a1),arrowsize);
		arrow->DrawArrow(x2-xfactor*aL*sin(a2),y2-aL*cos(a2),x2+xfactor*aL*sin(a2),y2+aL*cos(a2),arrowsize);
		arrow->DrawArrow(x3-xfactor*aL*sin(a3),y3-aL*cos(a3),x3+xfactor*aL*sin(a3),y3+aL*cos(a3),arrowsize);
		arrow->DrawArrow(x4-xfactor*aL*sin(a4),y4-aL*cos(a4),x4+xfactor*aL*sin(a4),y4+aL*cos(a4),arrowsize);
		arrow->DrawArrow(x5-xfactor*aL*sin(a5),y5-aL*cos(a5),x5+xfactor*aL*sin(a5),y5+aL*cos(a5),arrowsize);
		arrow->DrawArrow(x6-xfactor*aL*sin(a6),y6-aL*cos(a6),x6+xfactor*aL*sin(a6),y6+aL*cos(a6),arrowsize);
		arrow->DrawArrow(x7-xfactor*aL*sin(a7),y7-aL*cos(a7),x7+xfactor*aL*sin(a7),y7+aL*cos(a7),arrowsize);
		lineA->SetLineWidth(wA);
		lineA2->SetLineWidth(wA);
		lineB->SetLineWidth(wB);
		lineB2->SetLineWidth(wB);
		lineA->Draw("SAME");
		lineA2->Draw("SAME");
		lineB->Draw("SAME");
		lineB2->Draw("SAME");
		el1->Draw("SAME");
		el2->Draw("SAME");
		l1->Draw("SAME");
		l2->Draw("SAME");
		text_l1->DrawLatex(el1x,el1y+0.1,Form("HWPin: %.1f#circ",lambda1/degtorad));
		text_l2->DrawLatex(el2x,el2y+0.1,Form("HWPout: %.1f#circ",lambda2/degtorad));
	
		TPad* padB = new TPad("padB","",0,0.3,0.1,0.7);
		padB->Draw();
		padB->cd();
		gPad->SetLeftMargin(0.5);
		h1_B->Draw("HIST F");


		can->cd();
		TPad* pad1 = new TPad("pad1","",0.75,0,1.0,0.4);
		pad1->Draw();
		pad1->cd();
		g_AB_kick->GetXaxis()->SetLimits(-1,4);
		g_AB_kick->GetXaxis()->SetRangeUser(-1,4);
		g_AB_kick->GetYaxis()->SetRangeUser(-0.4,1.1);
		gPad->SetBottomMargin(0.2);
		gPad->SetLeftMargin(0.2);
		g_AB_kick->Draw("APL");
		gPad->SetGridy();

		can->Update();
		can->SaveAs(Form("gif_kick/%04d.png",framecount++));
		
	}

}
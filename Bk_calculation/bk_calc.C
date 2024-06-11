void bk_calc(){

	TFile* f1_R0 = TFile::Open("INFN_golden_R0_Bon.root");
	TGraph* g_trace_R0 = (TGraph*)f1_R0->Get("g_R0_Bon");

	TFile* f1_R1 = TFile::Open("INFN_golden_R1_Bon.root");
	TGraph* g_trace_R1 = (TGraph*)f1_R1->Get("g_R1_Bon");

	TFile* f2 = TFile::Open("Run4_S3010_B_corrected_eneBinned_FitEbinned_binAna_AsymWeighted.root");
	TH1D* h1_wiggle = (TH1D*)f2->Get("wiggle_plot_E1020");

	ifstream f_run3b;
	f_run3b.open("Beam_distribution_run3b.csv");
	f_run3b.ignore(256,'\n');
	TH2D* h2_run3b = new TH2D("h2_run3b","Beam distribution Run3b;x [mm];y [mm]",180,-45,45,180,-45,45);
	for (int i=0; i<180; i++){
		for (int j=0; j<180; j++){
			char comma;
			double x, y, z;
			f_run3b>>x>>comma>>y>>comma>>z;
			h2_run3b->Fill(x,y,z);
		}
	}
	f_run3b.close();
	TH1D* h1_run3b = h2_run3b->ProjectionX("h1_run3b");
	h1_run3b->Scale(1./h1_run3b->Integral());



	double fit_start = 30.1384; // us
	double fit_end = 650.064; // us

	//Zero the wiggle outside range
	for (int bx=0; bx<=h1_wiggle->GetNbinsX()+1; bx++){
		if (h1_wiggle->GetBinCenter(bx) < fit_start || h1_wiggle->GetBinCenter(bx) > fit_end){
			h1_wiggle->SetBinContent(bx,0);
		} else {
			//h1_wiggle->SetBinContent(bx,1);
		}
	}

	//Normalize wiggle
	double wiggle_integral = h1_wiggle->Integral();
	h1_wiggle->Scale(1./wiggle_integral);


	//Transform transient in ppm and in us
	double B = 1.45e7; // mG
	double f_azimuth = 0.085;
	double f_kickers = (53.1+53.0+55.0)/(3*55.0);


	//Transform transient in interpolated th1
	TH1D* h1_transient_Emma = new TH1D("h1_transient_Emma","Transient",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_transient_Emma->GetNbinsX(); bx++){
		double x = 0.001*h1_transient_Emma->GetBinCenter(bx);
		h1_transient_Emma->SetBinContent(bx,-35*exp(-x/0.0474)*f_azimuth*f_kickers*1e9/B);
		//h1_transient_Emma->SetBinContent(bx,-100);
		//h1_transient_Emma->SetBinContent(bx,x<0.05?-100:0);
	}

	TH1D* h1_transient_R0 = new TH1D("h1_transient_R0","Transient",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_transient_R0->GetNbinsX(); bx++){
		double x = 0.001*h1_transient_R0->GetBinCenter(bx);
		h1_transient_R0->SetBinContent(bx,(f_kickers*f_azimuth*1e9/B)*g_trace_R0->Eval(x));
	}

	TH1D* h1_transient_R1 = new TH1D("h1_transient_R1","Transient",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_transient_R1->GetNbinsX(); bx++){
		double x = 0.001*h1_transient_R1->GetBinCenter(bx);
		h1_transient_R1->SetBinContent(bx,(f_kickers*f_azimuth*1e9/B)*g_trace_R1->Eval(x));
	}


	TH1D* h1_runningavg_Emma = new TH1D("h1_runningavg_Emma","w_{a} bias = #frac{1}{t-30}#int_{30}^{t}Bk(t')dt'",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_runningavg_Emma->GetNbinsX(); bx++){
		int bin30 = h1_runningavg_Emma->FindBin(30);
		double binW = h1_runningavg_Emma->GetBinWidth(1);
		double t = h1_runningavg_Emma->GetBinCenter(bx);
		if (bx < bin30){
			h1_runningavg_Emma->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_Emma->Integral(bin30,bx);
			h1_runningavg_Emma->SetBinContent(bx,integral/(bx+1-bin30));
		}
	}
	TH1D* h1_runningavg_R0 = new TH1D("h1_runningavg_R0","w_{a} bias = #frac{1}{t-30}#int_{30}^{t}Bk(t')dt'",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_runningavg_R0->GetNbinsX(); bx++){
		int bin30 = h1_runningavg_R0->FindBin(30);
		double binW = h1_runningavg_R0->GetBinWidth(1);
		double t = h1_runningavg_R0->GetBinCenter(bx);
		if (bx < bin30){
			h1_runningavg_R0->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R0->Integral(bin30,bx);
			h1_runningavg_R0->SetBinContent(bx,integral/(bx+1-bin30));
		}
	}
	TH1D* h1_runningavg_R1 = new TH1D("h1_runningavg_R1","w_{a} bias = #frac{1}{t-30}#int_{30}^{t}Bk(t')dt'",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_runningavg_R1->GetNbinsX(); bx++){
		int bin30 = h1_runningavg_R1->FindBin(30);
		double binW = h1_runningavg_R1->GetBinWidth(1);
		double t = h1_runningavg_R1->GetBinCenter(bx);
		if (bx < bin30){
			h1_runningavg_R1->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R1->Integral(bin30,bx);
			h1_runningavg_R1->SetBinContent(bx,integral/(bx+1-bin30));
		}
	}


	TH1D* h1_A_R0 = new TH1D("h1_A_R0","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_A_R0->GetNbinsX(); bx++){
		int bin30 = h1_A_R0->FindBin(30);
		double binW = h1_A_R0->GetBinWidth(1);
		double t = h1_A_R0->GetBinCenter(bx);
		if (bx < bin30){
			h1_A_R0->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R0->Integral(1,bx,"width");
			h1_A_R0->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4)*(t-30));
		}
	}
	TH1D* h1_A_R1 = new TH1D("h1_A_R1","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_A_R1->GetNbinsX(); bx++){
		int bin30 = h1_A_R1->FindBin(30);
		double binW = h1_A_R1->GetBinWidth(1);
		double t = h1_A_R1->GetBinCenter(bx);
		if (bx < bin30){
			h1_A_R1->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R1->Integral(1,bx,"width");
			h1_A_R1->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4)*(t-30));
		}
	}
	TH1D* h1_A_Emma = new TH1D("h1_A_Emma","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_A_Emma->GetNbinsX(); bx++){
		int bin30 = h1_A_Emma->FindBin(30);
		double binW = h1_A_Emma->GetBinWidth(1);
		double t = h1_A_Emma->GetBinCenter(bx);
		if (bx < bin30){
			h1_A_Emma->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_Emma->Integral(1,bx,"width");
			h1_A_Emma->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4)*(t-30));
		}
	}

	TH1D* h1_B_R0 = new TH1D("h1_B_R0","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_B_R0->GetNbinsX(); bx++){
		int bin30 = h1_B_R0->FindBin(30);
		double binW = h1_B_R0->GetBinWidth(1);
		double t = h1_B_R0->GetBinCenter(bx);
		if (bx < bin30){
			h1_B_R0->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R0->Integral(1,bx,"width");
			h1_B_R0->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4));
		}
	}
	TH1D* h1_B_R1 = new TH1D("h1_B_R1","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_B_R1->GetNbinsX(); bx++){
		int bin30 = h1_B_R1->FindBin(30);
		double binW = h1_B_R1->GetBinWidth(1);
		double t = h1_B_R1->GetBinCenter(bx);
		if (bx < bin30){
			h1_B_R1->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_R1->Integral(1,bx,"width");
			h1_B_R1->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4));
		}
	}
	TH1D* h1_B_Emma = new TH1D("h1_B_Emma","A",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_B_Emma->GetNbinsX(); bx++){
		int bin30 = h1_B_Emma->FindBin(30);
		double binW = h1_B_Emma->GetBinWidth(1);
		double t = h1_B_Emma->GetBinCenter(bx);
		if (bx < bin30){
			h1_B_Emma->SetBinContent(bx,0);
		} else {
			double integral = h1_transient_Emma->Integral(1,bx,"width");
			h1_B_Emma->SetBinContent(bx,integral*pow(M_E,-(t-30)/64.4));
		}
	}

	//h1_wiggle->Rebin(30);
	//h1_transient->Rebin(30);
	//h1_transient->Scale(1./30);

	TH1D* h1_convolution_Emma = new TH1D("h1_convolution_Emma","B * wiggle",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_convolution_Emma->GetNbinsX(); bx++){
		h1_convolution_Emma->SetBinContent(bx,h1_runningavg_Emma->GetBinContent(bx)*h1_wiggle->GetBinContent(bx));
	}
	TH1D* h1_convolution_Emma_cumulative = (TH1D*)h1_convolution_Emma->GetCumulative();

	TH1D* h1_convolution_R0 = new TH1D("h1_convolution_R0","B * wiggle",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_convolution_R0->GetNbinsX(); bx++){
		h1_convolution_R0->SetBinContent(bx,h1_runningavg_R0->GetBinContent(bx)*h1_wiggle->GetBinContent(bx));
	}
	TH1D* h1_convolution_R0_cumulative = (TH1D*)h1_convolution_R0->GetCumulative();

	TH1D* h1_convolution_R1 = new TH1D("h1_convolution_R1","B * wiggle",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_convolution_R1->GetNbinsX(); bx++){
		h1_convolution_R1->SetBinContent(bx,h1_runningavg_R1->GetBinContent(bx)*h1_wiggle->GetBinContent(bx));
	}
	TH1D* h1_convolution_R1_cumulative = (TH1D*)h1_convolution_R1->GetCumulative();


	double xmin = 0;
	double xmax = 700;
	double ymin = -400;
	double ymax = 100;


	TLine* l0 = new TLine(xmin,0,xmax,0);
	l0->SetLineWidth(1);
	l0->SetLineColor(kBlack);

	TLine* l30 = new TLine(fit_start,ymin,fit_start,ymax);
	l30->SetLineWidth(2);
	l30->SetLineColor(kGreen);



	g_trace_R0->SetLineWidth(2);
	g_trace_R0->SetLineColor(kBlue);
	g_trace_R0->GetXaxis()->SetTitle("Time [#mus]");
	g_trace_R0->GetYaxis()->SetTitle("#Delta B [ppb]");
	g_trace_R1->SetLineWidth(2);
	g_trace_R1->SetLineColor(kRed);
	g_trace_R1->GetXaxis()->SetTitle("Time [#mus]");
	g_trace_R1->GetYaxis()->SetTitle("#Delta B [ppb]");

	h1_transient_Emma->SetLineWidth(2);
	h1_transient_Emma->SetLineColor(kBlack);
	h1_transient_Emma->GetXaxis()->SetTitle("Time [#mus]");
	h1_transient_Emma->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_transient_R0->SetLineWidth(2);
	h1_transient_R0->SetLineColor(kBlue);
	h1_transient_R0->GetXaxis()->SetTitle("Time [#mus]");
	h1_transient_R0->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_transient_R1->SetLineWidth(2);
	h1_transient_R1->SetLineColor(kRed);
	h1_transient_R1->GetXaxis()->SetTitle("Time [#mus]");
	h1_transient_R1->GetYaxis()->SetTitle("#Delta B [ppb]");

	h1_runningavg_Emma->SetLineWidth(2);
	h1_runningavg_Emma->SetLineColor(kBlack);
	h1_runningavg_Emma->GetXaxis()->SetTitle("Time [#mus]");
	h1_runningavg_Emma->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_runningavg_R0->SetLineWidth(2);
	h1_runningavg_R0->SetLineColor(kBlue);
	h1_runningavg_R0->GetXaxis()->SetTitle("Time [#mus]");
	h1_runningavg_R0->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_runningavg_R1->SetLineWidth(2);
	h1_runningavg_R1->SetLineColor(kRed);
	h1_runningavg_R1->GetXaxis()->SetTitle("Time [#mus]");
	h1_runningavg_R1->GetYaxis()->SetTitle("#Delta B [ppb]");

	h1_wiggle->SetTitle("Normalized wiggle [AMethod]");
	h1_wiggle->GetXaxis()->SetTitle("Time [#mus]");
	h1_wiggle->GetYaxis()->SetTitle("Positrons");

	h1_convolution_Emma->SetLineWidth(2);
	h1_convolution_Emma->SetLineColor(kBlack);
	h1_convolution_Emma->SetTitle("B \\circledast N");
	h1_convolution_Emma->GetXaxis()->SetTitle("Time [#mus]");
	h1_convolution_Emma->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_convolution_R0->SetLineWidth(2);
	h1_convolution_R0->SetLineColor(kBlue);
	h1_convolution_R0->SetTitle("B \\circledast N");
	h1_convolution_R0->GetXaxis()->SetTitle("Time [#mus]");
	h1_convolution_R0->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_convolution_R1->SetLineWidth(2);
	h1_convolution_R1->SetLineColor(kRed);
	h1_convolution_R1->SetTitle("B \\circledast N");
	h1_convolution_R1->GetXaxis()->SetTitle("Time [#mus]");
	h1_convolution_R1->GetYaxis()->SetTitle("#Delta B [ppb]");


	gStyle->SetOptStat(0);

	new TCanvas();
	g_trace_R0->Draw("AL");
	g_trace_R1->Draw("L");

	new TCanvas();
	h1_transient_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_transient_R0->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_transient_R0->Draw("HIST");
	h1_transient_R1->Draw("HIST SAME");
	h1_transient_Emma->Draw("HIST SAME");
	l0->Draw("SAME");
	l30->Draw("SAME");
	gPad->SetGridy();
	TLegend* leg0 = new TLegend(0.5,0.2,0.6,0.4);
	leg0->AddEntry(h1_transient_R0,"R0","L");
	leg0->AddEntry(h1_transient_R1,"R1","L");
	leg0->AddEntry(h1_transient_Emma,"Emma","L");
	leg0->Draw();

	TCanvas* can = new TCanvas("can","",1200,1200);
	can->Divide(1,4);

	can->cd(1);
	h1_transient_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_transient_R0->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_transient_R0->Draw("HIST");
	h1_transient_R1->Draw("HIST SAME");
	h1_transient_Emma->Draw("HIST SAME");
	l0->Draw("SAME");
	l30->Draw("SAME");
	gPad->SetGridy();
	TLegend* leg1 = new TLegend(0.5,0.2,0.6,0.4);
	leg1->AddEntry(h1_transient_R0,"R0","L");
	leg1->AddEntry(h1_transient_R1,"R1","L");
	leg1->Draw();
	cout<<"Transient R0 integral: "<<h1_transient_R0->Integral("width")<<"\n";
	cout<<"Transient R1 integral: "<<h1_transient_R1->Integral("width")<<"\n";


	can->cd(2);
	h1_runningavg_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_runningavg_R0->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_runningavg_R0->Draw("HIST");
	h1_runningavg_R1->Draw("HIST SAME");
	h1_runningavg_Emma->Draw("HIST SAME");
	l0->Draw("SAME");
	l30->Draw("SAME");
	gPad->SetGridy();
	TLegend* leg2 = new TLegend(0.5,0.2,0.6,0.4);
	leg2->AddEntry(h1_runningavg_R0,"R0","L");
	leg2->AddEntry(h1_runningavg_R1,"R1","L");
	leg2->Draw();
	cout<<"Runningavg R0 integral: "<<h1_runningavg_R0->Integral("width")<<"\n";
	cout<<"Runningavg R1 integral: "<<h1_runningavg_R1->Integral("width")<<"\n";


	can->cd(3);
	h1_wiggle->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_wiggle->Draw("HIST");
	gPad->SetGridy();
	cout<<"Wiggle integral: "<<h1_wiggle->Integral()<<"\n";

	can->cd(4);
	h1_convolution_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_convolution_R0->GetYaxis()->SetRangeUser(2e-4*ymin,2e-4*ymax);
	h1_convolution_R0->Draw("HIST");
	h1_convolution_R1->Draw("HIST SAME");
	h1_convolution_Emma->Draw("HIST SAME");
	gPad->SetGridy();
	TLegend* leg3 = new TLegend(0.5,0.2,0.7,0.4);

	double Bk_Emma = h1_convolution_Emma->Integral();
	double Bk_R0 = h1_convolution_R0->Integral();
	double Bk_R1 = h1_convolution_R1->Integral();

	leg3->AddEntry(h1_convolution_R0,Form("R0 : %.1f ppb",Bk_R0),"L");
	leg3->AddEntry(h1_convolution_R1,Form("R1 : %.1f ppb",Bk_R1),"L");
	leg3->AddEntry(h1_convolution_Emma,Form("Const : %.1f ppb",Bk_Emma),"L");
	leg3->Draw();
	cout<<"Convolution R0 integral: "<<Bk_R0<<"\n";
	cout<<"Convolution R1 integral: "<<Bk_R1<<"\n";
	cout<<"Convolution Emma integral: "<<Bk_Emma<<"\n";


	new TCanvas();
	h1_convolution_R0_cumulative->Draw("HIST");
	h1_convolution_R1_cumulative->Draw("HIST SAME");

	TLine* l_R0 = new TLine();
	TLine* l_R1 = new TLine();

	TGraphErrors* g_Bk_x = new TGraphErrors();
	g_Bk_x->SetPoint(0,0,Bk_R0);
	g_Bk_x->SetPointError(0,2,0.05*abs(Bk_R0));
	g_Bk_x->SetPoint(1,17.5,Bk_R1);
	g_Bk_x->SetPointError(1,2,0.10*abs(Bk_R1));

	new TCanvas();
	g_Bk_x->GetXaxis()->SetLimits(-50,50);
	g_Bk_x->GetXaxis()->SetRangeUser(-50,50);
	g_Bk_x->GetYaxis()->SetRangeUser(-400,0);
	g_Bk_x->GetXaxis()->SetTitle("x [mm]");
	g_Bk_x->GetYaxis()->SetTitle("Bk [ppb]");
	g_Bk_x->SetMarkerStyle(20);
	g_Bk_x->Draw("AP");

	TF1* f_model = new TF1("f_model","[0]+[2]*(x-[1])*(x-[1])",-45,45);
	f_model->SetParameters(Bk_R0,0,-0.5);
	f_model->FixParameter(1,0);
	g_Bk_x->Fit(f_model,"N");
	f_model->Draw("SAME");

	TH1D* h1_Bk_x_conv = new TH1D("h1_Bk_x_conv","Bk \\circledast N;x [mm];Bk [ppb]",h1_run3b->GetNbinsX(),h1_run3b->GetXaxis()->GetXmin(),h1_run3b->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_Bk_x_conv->GetNbinsX(); bx++){
		double x = h1_Bk_x_conv->GetBinCenter(bx);
		h1_Bk_x_conv->SetBinContent(bx,h1_run3b->GetBinContent(bx)*f_model->Eval(x));
	}

	new TCanvas();
	h1_run3b->Draw("HIST");
	l_R0->DrawLine(0,h1_run3b->GetMinimum(),0,h1_run3b->GetMaximum());
	l_R1->DrawLine(17.5,h1_run3b->GetMinimum(),17.5,h1_run3b->GetMaximum());

	new TCanvas();
	h1_Bk_x_conv->Draw("HIST");
	l_R0->DrawLine(0,h1_Bk_x_conv->GetMinimum(),0,h1_Bk_x_conv->GetMaximum());
	l_R1->DrawLine(17.5,h1_Bk_x_conv->GetMinimum(),17.5,h1_Bk_x_conv->GetMaximum());
	double Bk_Rmodel = h1_Bk_x_conv->Integral();
	TLegend* leg_Bkx = new TLegend(0.3,0.7,0.7,0.8);
	leg_Bkx->AddEntry(h1_Bk_x_conv,Form("Integral : %.1f ppb",Bk_Rmodel),"L");
	leg_Bkx->Draw();
	cout<<"Convolution space integral: "<<h1_Bk_x_conv->Integral()<<"\n";


	new TCanvas();
	h1_A_R0->SetLineColor(kBlue);
	h1_B_R0->SetLineColor(kRed);
	h1_A_R0->Draw("HIST");
	h1_B_R0->Draw("HIST SAME");
	new TCanvas();
	h1_A_R1->SetLineColor(kBlue);
	h1_B_R1->SetLineColor(kRed);
	h1_A_R1->Draw("HIST");
	h1_B_R1->Draw("HIST SAME");
	new TCanvas();
	h1_A_Emma->SetLineColor(kBlue);
	h1_B_Emma->SetLineColor(kRed);
	h1_A_Emma->Draw("HIST");
	h1_B_Emma->Draw("HIST SAME");

	cout<<"R0:\n";
	double integral_A_R0 = h1_A_R0->Integral(h1_A_R0->FindBin(0),h1_A_R0->FindBin(650),"width")/(64.4*64.4*64.4);	
	double integral_B_R0 = h1_B_R0->Integral(h1_B_R0->FindBin(0),h1_B_R0->FindBin(650),"width")/(64.4*64.4);
	cout<<"A integral = "<<integral_A_R0<<"\n";
	cout<<"B integral = "<<integral_B_R0<<"\n";
	double wa_bias_R0 = integral_A_R0 - integral_B_R0;
	cout<<"wa bias = "<<wa_bias_R0<<"\n";

	cout<<"R1:\n";
	double integral_A_R1 = h1_A_R1->Integral(h1_A_R1->FindBin(0),h1_A_R1->FindBin(650),"width")/(64.4*64.4*64.4);	
	double integral_B_R1 = h1_B_R1->Integral(h1_B_R1->FindBin(0),h1_B_R1->FindBin(650),"width")/(64.4*64.4);
	cout<<"A integral = "<<integral_A_R1<<"\n";
	cout<<"B integral = "<<integral_B_R1<<"\n";
	double wa_bias_R1 = integral_A_R1 - integral_B_R1;
	cout<<"wa bias = "<<wa_bias_R1<<"\n";

	cout<<"Emma:\n";
	double integral_A_Emma = h1_A_Emma->Integral(h1_A_Emma->FindBin(0),h1_A_Emma->FindBin(650),"width")/(64.4*64.4*64.4);	
	double integral_B_Emma = h1_B_Emma->Integral(h1_B_Emma->FindBin(0),h1_B_Emma->FindBin(650),"width")/(64.4*64.4);
	cout<<"A integral = "<<integral_A_Emma<<"\n";
	cout<<"B integral = "<<integral_B_Emma<<"\n";
	double wa_bias_Emma = integral_A_Emma - integral_B_Emma;
	cout<<"wa bias = "<<wa_bias_Emma<<"\n";

}
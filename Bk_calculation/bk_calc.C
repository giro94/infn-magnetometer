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


	//h1_wiggle->Rebin(30);
	//h1_transient->Rebin(30);
	//h1_transient->Scale(1./30);

	TH1D* h1_convolution_R0 = new TH1D("h1_convolution_R0","B * wiggle",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_convolution_R0->GetNbinsX(); bx++){
		h1_convolution_R0->SetBinContent(bx,h1_transient_R0->GetBinContent(bx)*h1_wiggle->GetBinContent(bx));
	}
	TH1D* h1_convolution_R0_cumulative = (TH1D*)h1_convolution_R0->GetCumulative();

	TH1D* h1_convolution_R1 = new TH1D("h1_convolution_R1","B * wiggle",h1_wiggle->GetNbinsX(),h1_wiggle->GetXaxis()->GetXmin(),h1_wiggle->GetXaxis()->GetXmax());
	for (int bx=1; bx<=h1_convolution_R1->GetNbinsX(); bx++){
		h1_convolution_R1->SetBinContent(bx,h1_transient_R1->GetBinContent(bx)*h1_wiggle->GetBinContent(bx));
	}
	TH1D* h1_convolution_R1_cumulative = (TH1D*)h1_convolution_R1->GetCumulative();


	double xmin = 0;
	double xmax = 300;
	double ymin = -600;
	double ymax = 100;

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

	h1_transient_R0->SetLineWidth(2);
	h1_transient_R0->SetLineColor(kBlue);
	h1_transient_R0->GetXaxis()->SetTitle("Time [#mus]");
	h1_transient_R0->GetYaxis()->SetTitle("#Delta B [ppb]");
	h1_transient_R1->SetLineWidth(2);
	h1_transient_R1->SetLineColor(kRed);
	h1_transient_R1->GetXaxis()->SetTitle("Time [#mus]");
	h1_transient_R1->GetYaxis()->SetTitle("#Delta B [ppb]");

	h1_wiggle->SetTitle("Normalized wiggle [AMethod]");
	h1_wiggle->GetXaxis()->SetTitle("Time [#mus]");
	h1_wiggle->GetYaxis()->SetTitle("Positrons");

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
	l30->Draw("SAME");

	TCanvas* can = new TCanvas("can","",1200,1200);
	can->Divide(1,3);

	can->cd(1);
	h1_transient_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_transient_R0->GetYaxis()->SetRangeUser(ymin,ymax);
	h1_transient_R0->Draw("HIST");
	h1_transient_R1->Draw("HIST SAME");
	l30->Draw("SAME");
	gPad->SetGridy();
	TLegend* leg1 = new TLegend(0.5,0.2,0.6,0.4);
	leg1->AddEntry(h1_transient_R0,"R0","L");
	leg1->AddEntry(h1_transient_R1,"R1","L");
	leg1->Draw();
	cout<<"Transient R0 integral: "<<h1_transient_R0->Integral()<<"\n";
	cout<<"Transient R1 integral: "<<h1_transient_R1->Integral()<<"\n";

	can->cd(2);
	h1_wiggle->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_wiggle->Draw("HIST");
	gPad->SetGridy();
	cout<<"Wiggle integral: "<<h1_wiggle->Integral()<<"\n";

	can->cd(3);
	h1_convolution_R0->GetXaxis()->SetRangeUser(xmin,xmax);
	h1_convolution_R0->GetYaxis()->SetRangeUser(1e-3*ymin,1e-3*ymax);
	h1_convolution_R0->Draw("HIST");
	h1_convolution_R1->Draw("HIST SAME");
	gPad->SetGridy();
	TLegend* leg3 = new TLegend(0.5,0.2,0.7,0.4);

	double Bk_R0 = h1_convolution_R0->Integral();
	double Bk_R1 = h1_convolution_R1->Integral();

	leg3->AddEntry(h1_convolution_R0,Form("R0 : %.1f ppb",Bk_R0),"L");
	leg3->AddEntry(h1_convolution_R1,Form("R1 : %.1f ppb",Bk_R1),"L");
	leg3->Draw();
	cout<<"Convolution R0 integral: "<<Bk_R0<<"\n";
	cout<<"Convolution R1 integral: "<<Bk_R1<<"\n";


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

	cout<<"Convolution space integral: "<<h1_Bk_x_conv->Integral()<<"\n";



}
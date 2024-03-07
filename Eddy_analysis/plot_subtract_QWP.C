#include "../analysis_tools.C"

void plot_subtract_QWP(){

	TFile* f1 = TFile::Open("analysis/analysis_EC_jan26_B5173_H25_Q130.root");
	TFile* f2 = TFile::Open("analysis/analysis_EC_jan28_B5173_H25_Q00.root");
	TFile* f3 = TFile::Open("analysis/analysis_EC_jan25_B5173_H25_Q22.5.root");


	TH1D* (*smoothing)(TH1D*,TString) = &runningAverage_5_10_15;
	double blum_norm_x = -0.323;
	double blum_norm_y = 129.;

	double blum_norm_x_trace = 4.769;

	TH1D* h1_kick1_p = ((TProfile*)f1->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_n = ((TProfile*)f2->Get("trace_kick1"))->ProjectionX();
	TH1D* h1_kick1_0 = ((TProfile*)f3->Get("trace_kick1"))->ProjectionX();

	TH1D* h1_kick1_p_calib = ((TProfile*)f1->Get("trace_kick1_calibrated"))->ProjectionX();
	TH1D* h1_kick1_n_calib = ((TProfile*)f2->Get("trace_kick1_calibrated"))->ProjectionX();
	h1_kick1_n_calib->Scale(-1);

	TH1D* h1_trace_p = ((TProfile*)f1->Get("trace"))->ProjectionX();
	TH1D* h1_trace_n = ((TProfile*)f2->Get("trace"))->ProjectionX();
	TH1D* h1_trace_0 = ((TProfile*)f3->Get("trace"))->ProjectionX();




	/////////////////////////////////////////////////


	TH1D* h1_kick1_pn_0 = (TH1D*)h1_kick1_p->Clone("h1_kick1_pn_0");
	h1_kick1_pn_0->Add(h1_kick1_n);
	h1_kick1_pn_0->Scale(0.5);
	h1_kick1_pn_0->Add(h1_kick1_0,-1);

	TH1D* h1_kick1_p_n = (TH1D*)h1_kick1_p->Clone("h1_kick1_p_n");
	h1_kick1_p_n->Add(h1_kick1_n,-1);

	TH1D* h1_kick1_n_p = (TH1D*)h1_kick1_n->Clone("h1_kick1_n_p");
	h1_kick1_n_p->Add(h1_kick1_p,-1);

	TH1D* h1_kick1_p_0 = (TH1D*)h1_kick1_p->Clone("h1_kick1_p_0");
	h1_kick1_p_0->Add(h1_kick1_0,-1);

	TH1D* h1_kick1_n_0 = (TH1D*)h1_kick1_n->Clone("h1_kick1_n_0");
	h1_kick1_n_0->Add(h1_kick1_0,-1);

	TH1D* h1_kick1_p0_n0 = (TH1D*)h1_kick1_p_0->Clone("h1_kick1_p0_n0");
	h1_kick1_p0_n0->Add(h1_kick1_n_0);
	h1_kick1_p0_n0->Scale(0.5);


	TH1D* h1_kick1_p_n_0 = (TH1D*)h1_kick1_p->Clone("h1_kick1_p_n_0");
	h1_kick1_p_n_0->Add(h1_kick1_n,-1);
	h1_kick1_p_n_0->Add(h1_kick1_0,0.5);


	TH1D* h1_kick1_pnc_0 = (TH1D*)h1_kick1_p_calib->Clone("h1_kick1_pnc_0");
	h1_kick1_pnc_0->Add(h1_kick1_n_calib);
	h1_kick1_pnc_0->Scale(0.5);
	h1_kick1_pnc_0->Add(h1_kick1_0,-1);

	TH1D* h1_kick1_pc_nc = (TH1D*)h1_kick1_p_calib->Clone("h1_kick1_pc_nc");
	h1_kick1_pc_nc->Add(h1_kick1_n_calib,-1);

	TH1D* h1_kick1_nc_pc = (TH1D*)h1_kick1_n_calib->Clone("h1_kick1_nc_pc");
	h1_kick1_nc_pc->Add(h1_kick1_p_calib,-1);

	TH1D* h1_kick1_pc_0 = (TH1D*)h1_kick1_p_calib->Clone("h1_kick1_pc_0");
	h1_kick1_pc_0->Add(h1_kick1_0,-1);

	TH1D* h1_kick1_nc_0 = (TH1D*)h1_kick1_n_calib->Clone("h1_kick1_nc_0");
	h1_kick1_nc_0->Add(h1_kick1_0,-1);



	TH1D* h1_kick1_pc_nc_0 = (TH1D*)h1_kick1_p_calib->Clone("h1_kick1_pc_nc_0");
	h1_kick1_pc_nc_0->Add(h1_kick1_n_calib,-1);
	h1_kick1_pc_nc_0->Scale(0.5);
	h1_kick1_pc_nc_0->Add(h1_kick1_0);




	TH1D** h1_kick1_pX_nY = new TH1D* [41];
	TH1D** h1_kick1_pX_nY_ra = new TH1D* [41];
	TGraph* g_rms = new TGraph();
	for (int i=0; i<41; i++){
		double x = i*0.025;
		double y = -(1-x);
		h1_kick1_pX_nY[i] = (TH1D*)h1_kick1_p_calib->Clone(Form("h1_kick1_p%.3f_n%.3f",x,y));
		h1_kick1_pX_nY[i]->Add(h1_kick1_p_calib,h1_kick1_n_calib,x,y);
		h1_kick1_pX_nY_ra[i] = smoothing(h1_kick1_pX_nY[i],"");
		h1_kick1_pX_nY_ra[i]->SetLineWidth(2);
		h1_kick1_pX_nY_ra[i]->GetXaxis()->SetRangeUser(2,6);
		h1_kick1_pX_nY_ra[i]->Scale(blum_norm_y/h1_kick1_pX_nY_ra[i]->Interpolate(blum_norm_x));
		TGraphErrors *g = new TGraphErrors(SetupFFT(h1_kick1_pX_nY_ra[i],2,6));
		double stddev = g->GetRMS(2);
		cout<<stddev<<"\n";
		g_rms->SetPoint(i,x,stddev);
	}



	TH1D** h1_kick1_pX_nY_k0 = new TH1D* [41];
	TH1D** h1_kick1_pX_nY_k0_ra = new TH1D* [41];
	TGraph* g_rms2 = new TGraph();
	for (int i=0; i<41; i++){
		double x = i*0.025;
		h1_kick1_pX_nY_k0[i] = (TH1D*)h1_kick1_pX_nY[26]->Clone(Form("h1_kick1_p0.65_n0.35_k0_%.3f",x));
		h1_kick1_pX_nY_k0[i]->Add(h1_kick1_0,x);
		h1_kick1_pX_nY_k0_ra[i] = smoothing(h1_kick1_pX_nY_k0[i],"");
		h1_kick1_pX_nY_k0_ra[i]->SetLineWidth(2);
		h1_kick1_pX_nY_k0_ra[i]->GetXaxis()->SetRangeUser(2,6);
		h1_kick1_pX_nY_k0_ra[i]->Scale(blum_norm_y/h1_kick1_pX_nY_k0_ra[i]->Interpolate(blum_norm_x));
		TGraphErrors *g2 = new TGraphErrors(SetupFFT(h1_kick1_pX_nY_k0_ra[i],2,6));
		double stddev = g2->GetRMS(2);
		cout<<stddev<<"\n";
		g_rms2->SetPoint(i,x,stddev);
	}


	TH1D*** h1_kick1_pX_nY_k0Z = new TH1D** [41];
	TH1D*** h1_kick1_pX_nY_k0Z_ra = new TH1D** [41];
	TH2D* g_rms3_ra = new TH2D("g_rms3_ra","",41,-0.0125,1.0125,41,-0.0125,1.0125);
	TGraph* g_rms3_ra_best = new TGraph();
	int best_i = 0;
	int best_j = 0;
	double best_x = 0;
	double best_y = 0;
	double best_z = 0;
	double best_rms3_ra = 1e6;
	for (int i=0; i<41; i++){
		double x = i*0.025;
		double y = -(1-x);
		h1_kick1_pX_nY_k0Z[i] = new TH1D* [41];
		h1_kick1_pX_nY_k0Z_ra[i] = new TH1D* [41];
		for (int j=0; j<41; j++){
			double z = j*0.025;
			h1_kick1_pX_nY_k0Z[i][j] = (TH1D*)h1_kick1_p->Clone(Form("h1_kick1_p%.3f_n%.3f_k0%.3f",x,y,z));
			h1_kick1_pX_nY_k0Z[i][j]->Add(h1_kick1_p,h1_kick1_n,x,y);
			h1_kick1_pX_nY_k0Z[i][j]->Add(h1_kick1_0,z);
			h1_kick1_pX_nY_k0Z_ra[i][j] = smoothing(h1_kick1_pX_nY_k0Z[i][j],"");
			h1_kick1_pX_nY_k0Z_ra[i][j]->SetLineWidth(2);
			h1_kick1_pX_nY_k0Z_ra[i][j]->GetXaxis()->SetRangeUser(2,6);
			h1_kick1_pX_nY_k0Z_ra[i][j]->Scale(blum_norm_y/h1_kick1_pX_nY_k0Z_ra[i][j]->Interpolate(blum_norm_x));

			TGraphErrors *g3_ra = new TGraphErrors(SetupFFT(h1_kick1_pX_nY_k0Z_ra[i][j],2,6));
			double stddev_ra = g3_ra->GetRMS(2);
			cout<<i<<" "<<j<<" "<<x<<" "<<y<<" "<<z<<" "<<stddev_ra<<"\n";
			g_rms3_ra->Fill(x,z,stddev_ra);
			if (stddev_ra < best_rms3_ra){
				best_rms3_ra = stddev_ra;
				best_i = i;
				best_j = j;
				best_x = x;
				best_y = y;
				best_z = z;
				g_rms3_ra_best->SetPoint(0,x,z);
			}
		}
	}


	TH1D* h1_trace_pX_nY_k0Z = (TH1D*)h1_trace_p->Clone(Form("h1_trace_p%.3f_n%.3f_k0%.3f",best_x,best_y,best_z));
	h1_trace_pX_nY_k0Z->Add(h1_trace_p,h1_trace_n,best_x,best_y);
	h1_trace_pX_nY_k0Z->Add(h1_trace_0,best_z);
	TH1D* h1_trace_pX_nY_k0Z_ra = smoothing(h1_trace_pX_nY_k0Z,"");
	h1_trace_pX_nY_k0Z_ra->SetLineWidth(2);
	h1_trace_pX_nY_k0Z_ra->Scale(blum_norm_y/h1_trace_pX_nY_k0Z_ra->Interpolate(blum_norm_x_trace));





	TH1D*** h1_kick1_pcX_ncY_k0Z = new TH1D** [41];
	TH1D*** h1_kick1_pcX_ncY_k0Z_ra = new TH1D** [41];
	TH2D* g_rms3_calib_ra = new TH2D("g_rms3_calib_ra","",41,-0.0125,1.0125,41,-0.0125,1.0125);
	TGraph* g_rms3_calib_ra_best = new TGraph();
	int best_calib_i = 0;
	int best_calib_j = 0;
	double best_calib_x = 0;
	double best_calib_z = 0;
	double best_rms3_calib_ra = 1e6;
	for (int i=0; i<41; i++){
		double x = i*0.025;
		double y = -(1-x);
		h1_kick1_pcX_ncY_k0Z[i] = new TH1D* [41];
		h1_kick1_pcX_ncY_k0Z_ra[i] = new TH1D* [41];
		for (int j=0; j<41; j++){
			double z = j*0.025;
			h1_kick1_pcX_ncY_k0Z[i][j] = (TH1D*)h1_kick1_p_calib->Clone(Form("h1_kick1_pc%.3f_nc%.3f_k0%.3f",x,y,z));
			h1_kick1_pcX_ncY_k0Z[i][j]->Add(h1_kick1_p_calib,h1_kick1_n_calib,x,y);
			h1_kick1_pcX_ncY_k0Z[i][j]->Add(h1_kick1_0,z);
			h1_kick1_pcX_ncY_k0Z_ra[i][j] = smoothing(h1_kick1_pcX_ncY_k0Z[i][j],"");
			h1_kick1_pcX_ncY_k0Z_ra[i][j]->SetLineWidth(2);
			h1_kick1_pcX_ncY_k0Z_ra[i][j]->GetXaxis()->SetRangeUser(2,6);
			h1_kick1_pcX_ncY_k0Z_ra[i][j]->Scale(blum_norm_y/h1_kick1_pcX_ncY_k0Z_ra[i][j]->Interpolate(blum_norm_x));


			TGraphErrors *g3_ra = new TGraphErrors(SetupFFT(h1_kick1_pcX_ncY_k0Z_ra[i][j],2,6));
			double stddev_ra = g3_ra->GetRMS(2);
			cout<<i<<" "<<j<<" "<<x<<" "<<y<<" "<<z<<" "<<stddev_ra<<"\n";
			g_rms3_calib_ra->Fill(x,z,stddev_ra);
			if (stddev_ra < best_rms3_calib_ra){
				best_rms3_calib_ra = stddev_ra;
				best_calib_i = i;
				best_calib_j = j;
				best_calib_x = x;
				best_calib_z = z;
				g_rms3_calib_ra_best->SetPoint(0,x,z);
			}
		}
	}



	///////////////////////// Smoothing



	TH1D* h1_kick1_p_ra = smoothing(h1_kick1_p,"");
	TH1D* h1_kick1_n_ra = smoothing(h1_kick1_n,"");
	TH1D* h1_kick1_0_ra = smoothing(h1_kick1_0,"");
	TH1D* h1_kick1_p_calib_ra = smoothing(h1_kick1_p_calib,"");
	TH1D* h1_kick1_n_calib_ra = smoothing(h1_kick1_n_calib,"");


	TH1D* h1_kick1_pn_0_ra = smoothing(h1_kick1_pn_0,"");
	TH1D* h1_kick1_p_n_ra = smoothing(h1_kick1_p_n,"");
	TH1D* h1_kick1_n_p_ra = smoothing(h1_kick1_n_p,"");
	TH1D* h1_kick1_p_0_ra = smoothing(h1_kick1_p_0,"");
	TH1D* h1_kick1_n_0_ra = smoothing(h1_kick1_n_0,"");
	TH1D* h1_kick1_p0_n0_ra = smoothing(h1_kick1_p0_n0,"");

	TH1D* h1_kick1_p_n_0_ra = smoothing(h1_kick1_p_n_0,"");


	TH1D* h1_kick1_pnc_0_ra = smoothing(h1_kick1_pnc_0,"");
	TH1D* h1_kick1_pc_nc_ra = smoothing(h1_kick1_pc_nc,"");
	TH1D* h1_kick1_nc_pc_ra = smoothing(h1_kick1_nc_pc,"");
	TH1D* h1_kick1_pc_0_ra = smoothing(h1_kick1_pc_0,"");
	TH1D* h1_kick1_nc_0_ra = smoothing(h1_kick1_nc_0,"");

	TH1D* h1_kick1_pc_nc_0_ra = smoothing(h1_kick1_pc_nc_0,"");


	///////////////////////// Normalizing


	h1_kick1_pn_0_ra->Scale(blum_norm_y/h1_kick1_pn_0_ra->Interpolate(blum_norm_x));
	h1_kick1_p_n_ra->Scale(blum_norm_y/h1_kick1_p_n_ra->Interpolate(blum_norm_x));
	h1_kick1_n_p_ra->Scale(blum_norm_y/h1_kick1_n_p_ra->Interpolate(blum_norm_x));
	h1_kick1_p_0_ra->Scale(blum_norm_y/h1_kick1_p_0_ra->Interpolate(blum_norm_x));
	h1_kick1_n_0_ra->Scale(blum_norm_y/h1_kick1_n_0_ra->Interpolate(blum_norm_x));
	h1_kick1_p0_n0_ra->Scale(blum_norm_y/h1_kick1_p0_n0_ra->Interpolate(blum_norm_x));

	h1_kick1_p_n_0_ra->Scale(blum_norm_y/h1_kick1_p_n_0_ra->Interpolate(blum_norm_x));

	h1_kick1_pnc_0_ra->Scale(blum_norm_y/h1_kick1_pnc_0_ra->Interpolate(blum_norm_x));
	h1_kick1_pc_nc_ra->Scale(blum_norm_y/h1_kick1_pc_nc_ra->Interpolate(blum_norm_x));
	h1_kick1_nc_pc_ra->Scale(blum_norm_y/h1_kick1_nc_pc_ra->Interpolate(blum_norm_x));
	h1_kick1_pc_0_ra->Scale(blum_norm_y/h1_kick1_pc_0_ra->Interpolate(blum_norm_x));
	h1_kick1_nc_0_ra->Scale(blum_norm_y/h1_kick1_nc_0_ra->Interpolate(blum_norm_x));

	h1_kick1_pc_nc_0_ra->Scale(blum_norm_y/h1_kick1_pc_nc_0_ra->Interpolate(blum_norm_x));



	//////////////////////////////// Drawing

	gStyle->SetPalette(kRainBow);
	gStyle->SetOptStat(0);



	h1_kick1_p_ra->SetLineColor(kBlack);
	h1_kick1_n_ra->SetLineColor(kBlue);
	h1_kick1_0_ra->SetLineColor(kRed);
	h1_kick1_p_calib_ra->SetLineColor(kBlack);
	h1_kick1_n_calib_ra->SetLineColor(kBlue);
	h1_kick1_p_ra->SetLineWidth(2);
	h1_kick1_n_ra->SetLineWidth(2);
	h1_kick1_0_ra->SetLineWidth(2);
	h1_kick1_p_calib_ra->SetLineWidth(2);
	h1_kick1_n_calib_ra->SetLineWidth(2);
	h1_kick1_p_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_n_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p_calib_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_n_calib_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_n_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p_calib_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_n_calib_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p_ra->SetTitle("QWP 130");
	h1_kick1_n_ra->SetTitle("QWP 0");
	h1_kick1_0_ra->SetTitle("QWP 22.5");
	h1_kick1_p_calib_ra->SetTitle("QWP 130");
	h1_kick1_n_calib_ra->SetTitle("QWP 0");


	h1_kick1_pn_0_ra->SetLineColor(1);
	h1_kick1_pn_0_ra->SetLineWidth(2);
	h1_kick1_pn_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_pn_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_pn_0_ra->SetTitle("(p+n)/2 - k0");

	h1_kick1_p_n_ra->SetLineColor(2);
	h1_kick1_p_n_ra->SetLineWidth(2);
	h1_kick1_p_n_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p_n_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p_n_ra->SetTitle("p-n");

	h1_kick1_n_p_ra->SetLineColor(2);
	h1_kick1_n_p_ra->SetLineWidth(2);
	h1_kick1_n_p_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_n_p_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_n_p_ra->SetTitle("n-p");

	h1_kick1_p_0_ra->SetLineColor(3);
	h1_kick1_p_0_ra->SetLineWidth(2);
	h1_kick1_p_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p_0_ra->SetTitle("p-k0");

	h1_kick1_n_0_ra->SetLineColor(4);
	h1_kick1_n_0_ra->SetLineWidth(2);
	h1_kick1_n_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_n_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_n_0_ra->SetTitle("n-k0");


	h1_kick1_p0_n0_ra->SetLineColor(6);
	h1_kick1_p0_n0_ra->SetLineWidth(2);
	h1_kick1_p0_n0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p0_n0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p0_n0_ra->SetTitle("0.4*(p-k0)-0.6*(n-k0)");


	h1_kick1_p_n_0_ra->SetLineColor(6);
	h1_kick1_p_n_0_ra->SetLineWidth(2);
	h1_kick1_p_n_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_p_n_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_p_n_0_ra->SetTitle("p-n+k0/2");


	h1_kick1_pnc_0_ra->SetLineColor(1);
	h1_kick1_pnc_0_ra->SetLineWidth(2);
	h1_kick1_pnc_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_pnc_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_pnc_0_ra->SetTitle("(pc+nc)/2 - k0");

	h1_kick1_pc_nc_ra->SetLineColor(2);
	h1_kick1_pc_nc_ra->SetLineWidth(2);
	h1_kick1_pc_nc_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_pc_nc_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_pc_nc_ra->SetTitle("pc-nc");

	h1_kick1_nc_pc_ra->SetLineColor(2);
	h1_kick1_nc_pc_ra->SetLineWidth(2);
	h1_kick1_nc_pc_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_nc_pc_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_nc_pc_ra->SetTitle("nc-pc");

	h1_kick1_pc_0_ra->SetLineColor(3);
	h1_kick1_pc_0_ra->SetLineWidth(2);
	h1_kick1_pc_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_pc_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_pc_0_ra->SetTitle("pc-k0");

	h1_kick1_nc_0_ra->SetLineColor(4);
	h1_kick1_nc_0_ra->SetLineWidth(2);
	h1_kick1_nc_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_nc_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_nc_0_ra->SetTitle("nc-k0");


	h1_kick1_pc_nc_0_ra->SetLineColor(6);
	h1_kick1_pc_nc_0_ra->SetLineWidth(2);
	h1_kick1_pc_nc_0_ra->GetXaxis()->SetRangeUser(2,4);
	h1_kick1_pc_nc_0_ra->GetYaxis()->SetRangeUser(-60,60);
	h1_kick1_pc_nc_0_ra->SetTitle("(pc-nc)/2+k0");



	new TCanvas();
	h1_kick1_p_ra->DrawCopy("HIST");
	h1_kick1_n_ra->DrawCopy("HIST SAME");
	h1_kick1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_p_calib_ra->DrawCopy("HIST");
	h1_kick1_n_calib_ra->DrawCopy("HIST SAME");
	h1_kick1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas();
	h1_kick1_pn_0_ra->DrawCopy("HIST");
	h1_kick1_p_n_ra->DrawCopy("HIST SAME");
	h1_kick1_p_0_ra->DrawCopy("HIST SAME");
	h1_kick1_n_0_ra->DrawCopy("HIST SAME");
	h1_kick1_p_n_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();



	new TCanvas();
	h1_kick1_pnc_0_ra->DrawCopy("HIST");
	h1_kick1_pc_nc_ra->DrawCopy("HIST SAME");
	h1_kick1_pc_0_ra->DrawCopy("HIST SAME");
	h1_kick1_nc_0_ra->DrawCopy("HIST SAME");
	h1_kick1_pc_nc_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();




	new TCanvas();
	h1_kick1_p_n_ra->DrawCopy("HIST");
	h1_kick1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_n_p_ra->DrawCopy("HIST");
	h1_kick1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_nc_pc_ra->Scale(0.5);
	h1_kick1_nc_pc_ra->DrawCopy("HIST");
	h1_kick1_0_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_p_0_ra->DrawCopy("HIST");
	h1_kick1_n_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	h1_kick1_n_0_ra->DrawCopy("HIST");
	h1_kick1_p_ra->DrawCopy("HIST SAME");
	gPad->BuildLegend();


	new TCanvas();
	for (int i=0; i<41; i++){
		h1_kick1_pX_nY_ra[i]->DrawCopy("HIST SAME PLC");
	}

	new TCanvas();
	g_rms->SetMarkerStyle(20);
	g_rms->Draw("APL");


	new TCanvas();
	h1_kick1_pX_nY_ra[26]->DrawCopy("HIST");
	h1_kick1_0_ra->DrawCopy("HIST SAME");


	new TCanvas();
	for (int i=0; i<41; i++){
		h1_kick1_pX_nY_k0_ra[i]->DrawCopy("HIST SAME PLC");
	}

	new TCanvas();
	g_rms2->SetMarkerStyle(20);
	g_rms2->Draw("APL");


	new TCanvas();
	h1_kick1_pX_nY_k0_ra[12]->DrawCopy("HIST");


	new TCanvas("","",1000,800);
	g_rms3_ra->SetTitle("Uncalibrated scan");
	g_rms3_ra->GetXaxis()->SetTitle("wPos #equiv 1-wNeg)");
	g_rms3_ra->GetYaxis()->SetTitle("wZero");
	g_rms3_ra->GetZaxis()->SetTitle("Trace RMS in (2,6) ms [mG]");
	g_rms3_ra->SetContour(100);
	g_rms3_ra->Draw("colz");
	g_rms3_ra_best->SetMarkerStyle(20);
	g_rms3_ra_best->SetMarkerColor(kRed);
	g_rms3_ra_best->Draw("P");
	TText* text_ra_best = new TText(best_x+0.01,best_z+0.01,Form("(%.2f,%.2f)",best_x,best_z));
	text_ra_best->SetTextSize(0.035);
	text_ra_best->SetTextColor(kRed);
	text_ra_best->Draw();


	new TCanvas("","",1000,800);
	g_rms3_calib_ra->SetTitle("Calibrated scan");
	g_rms3_calib_ra->GetXaxis()->SetTitle("wPos #equiv 1-wNeg)");
	g_rms3_calib_ra->GetYaxis()->SetTitle("wZero");
	g_rms3_calib_ra->GetZaxis()->SetTitle("Trace RMS in (2,6) ms [mG]");
	g_rms3_calib_ra->SetContour(100);
	g_rms3_calib_ra->Draw("colz");
	g_rms3_calib_ra_best->SetMarkerStyle(20);
	g_rms3_calib_ra_best->SetMarkerColor(kRed);
	g_rms3_calib_ra_best->Draw("P");
	TText* text_calib_ra_best = new TText(best_calib_x+0.01,best_calib_z+0.01,Form("(%.2f,%.2f)",best_calib_x,best_calib_z));
	text_calib_ra_best->SetTextSize(0.035);
	text_calib_ra_best->SetTextColor(kRed);
	text_calib_ra_best->Draw();

	new TCanvas();
	h1_kick1_p_calib_ra->GetXaxis()->SetRangeUser(-0.6,2.0);
	h1_kick1_p_calib_ra->GetYaxis()->SetRangeUser(-50,150);
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->GetXaxis()->SetRangeUser(-0.6,2.0);
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->GetYaxis()->SetRangeUser(-50,150);
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->SetLineColor(6);
	h1_kick1_p_n_0_ra->GetXaxis()->SetRangeUser(-0.6,2.0);
	h1_kick1_p_n_0_ra->GetYaxis()->SetRangeUser(-50,150);
	h1_kick1_p_calib_ra->DrawCopy("HIST");
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->DrawCopy("HIST SAME");
	//h1_kick1_p_n_0_ra->DrawCopy("HIST SAME");
	TLegend* leg = new TLegend(0.5,0.6,0.85,0.8);
	leg->AddEntry(h1_kick1_p_calib_ra,"Pos","L");
	leg->AddEntry(h1_kick1_pX_nY_k0Z_ra[best_i][best_j],"0.6*Pos -0.4*Neg +0.425*Zero","L");
	//leg->AddEntry(h1_kick1_p_n_0_ra,"0.54*Pos -0.46*Neg +0.5*Zero","L");
	leg->Draw();
	gPad->SetGridy();


	new TCanvas();
	h1_kick1_pc_nc_0_ra->SetLineColor(kBlue);
	h1_kick1_pc_nc_0_ra->DrawCopy("HIST");
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_p_calib->DrawCopy("HIST");
	h1_kick1_pX_nY_k0Z[best_i][best_j]->DrawCopy("HIST SAME");
	h1_kick1_p_n_0->DrawCopy("HIST SAME");



	new TCanvas();
	h1_kick1_p_0_ra->DrawCopy("HIST");
	h1_kick1_n_0_ra->DrawCopy("HIST SAME");
	h1_kick1_p0_n0_ra->DrawCopy("HIST SAME");
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->DrawCopy("HIST SAME");
	gPad->BuildLegend();

	new TCanvas();
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->DrawCopy("HIST");
	h1_kick1_pcX_ncY_k0Z_ra[best_calib_i][best_calib_j]->DrawCopy("HIST SAME");


	new TCanvas();
	h1_trace_pX_nY_k0Z_ra->Draw("HIST");


	new TCanvas();
	h1_kick1_p_0_ra->DrawCopy("HIST");
	h1_kick1_p_n_ra->DrawCopy("HIST SAME");



	new TCanvas();
	h1_kick1_pn_0_ra->SetLineColor(kBlack);
	h1_kick1_pn_0_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_pn_0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_pn_0_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_pn_0_ra->DrawCopy("HIST");
	gPad->SetGridy();

	new TCanvas();
	h1_kick1_p_n_ra->SetLineColor(kBlack);
	h1_kick1_p_n_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_p_n_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_p_n_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_p_n_ra->DrawCopy("HIST");
	gPad->SetGridy();


	new TCanvas();
	h1_kick1_p_0_ra->SetLineColor(kBlack);
	h1_kick1_p_0_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_p_0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_p_0_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_p_0_ra->DrawCopy("HIST");
	gPad->SetGridy();

	new TCanvas();
	h1_kick1_n_0_ra->SetLineColor(kBlack);
	h1_kick1_n_0_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_n_0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_n_0_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_n_0_ra->DrawCopy("HIST");
	gPad->SetGridy();

	new TCanvas();
	h1_kick1_p_n_0_ra->SetLineColor(kBlack);
	h1_kick1_p_n_0_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_p_n_0_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_p_n_0_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_p_n_0_ra->DrawCopy("HIST");
	gPad->SetGridy();

	new TCanvas();
	h1_kick1_p_calib_ra->SetLineColor(kBlack);
	h1_kick1_p_calib_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_p_calib_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_p_calib_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_p_calib_ra->DrawCopy("HIST");
	gPad->SetGridy();

	new TCanvas();
	h1_kick1_n_calib_ra->Scale(-1);
	h1_kick1_n_calib_ra->SetLineColor(kBlack);
	h1_kick1_n_calib_ra->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_n_calib_ra->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_n_calib_ra->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_n_calib_ra->DrawCopy("HIST");
	gPad->SetGridy();


	new TCanvas();
	h1_kick1_pX_nY_k0Z[best_i][best_j]->SetLineColor(kBlack);
	h1_kick1_pX_nY_k0Z[best_i][best_j]->GetXaxis()->SetRangeUser(-0.6,3.0);
	h1_kick1_pX_nY_k0Z[best_i][best_j]->GetYaxis()->SetRangeUser(-50,50);
	h1_kick1_pX_nY_k0Z[best_i][best_j]->GetYaxis()->SetTitle("Bfield [mG]");
	h1_kick1_pX_nY_k0Z[best_i][best_j]->DrawCopy("HIST");
	h1_kick1_pX_nY_k0Z_ra[best_i][best_j]->DrawCopy("HIST SAME");
	gPad->SetGridy();
}
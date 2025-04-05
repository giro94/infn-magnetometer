void convert_UMass_root(){

	ifstream f0,f1,f2,f3,f4,f5,f6,f7;

	f0.open("AveragedKicks_Times.csv");
	f1.open("AveragedKicks_Cali_fullFieldSec1.csv");
	f2.open("AveragedKicks_Cali_fullFieldSec2.csv");
	f3.open("AveragedKicks_Cali_fullFieldSec3.csv");
	f4.open("AveragedKicks_Cali_fullFieldSec4.csv");
	f5.open("AveragedKicks_Cali_fullFieldSec5.csv");


	f6.open("AveragedKicksTimes_Total2022.csv");
	f7.open("MaxAveragedKicks_Total2022.csv");

	TFile* fout = new TFile("UMass.root","recreate");

	TGraph* g_K1_R0 = new TGraph();
	TGraph* g_K1_R3p2 = new TGraph();
	TGraph* g_K1_R6p6 = new TGraph();
	TGraph* g_K1_Rm6p6 = new TGraph();
	TGraph* g_K3_R0 = new TGraph();
	TGraph* g_K1_R0_2022 = new TGraph();
	
	g_K1_R0->SetName("g_K1_R0");
	g_K1_R3p2->SetName("g_K1_R3p2");
	g_K1_R6p6->SetName("g_K1_R6p6");
	g_K1_Rm6p6->SetName("g_K1_Rm6p6");
	g_K3_R0->SetName("g_K3_R0");
	g_K1_R0_2022->SetName("g_K1_R0_2022");

	g_K1_R0->SetTitle("K1 0.0 mm");
	g_K1_R3p2->SetTitle("K1 3.2 mm");
	g_K1_R6p6->SetTitle("K1 6.6 mm");
	g_K1_Rm6p6->SetTitle("K1 -6.6 mm");
	g_K3_R0->SetTitle("K3 0.0 mm");
	g_K1_R0_2022->SetTitle("K1 2022 0.0 mm");

	g_K1_R0->GetXaxis()->SetTitle("Time [ms]");
	g_K1_R3p2->GetXaxis()->SetTitle("Time [ms]");
	g_K1_R6p6->GetXaxis()->SetTitle("Time [ms]");
	g_K1_Rm6p6->GetXaxis()->SetTitle("Time [ms]");
	g_K3_R0->GetXaxis()->SetTitle("Time [ms]");
	g_K1_R0_2022->GetXaxis()->SetTitle("Time [ms]");

	g_K1_R0->GetYaxis()->SetTitle("Trace [mG]");
	g_K1_R3p2->GetYaxis()->SetTitle("Trace [mG]");
	g_K1_R6p6->GetYaxis()->SetTitle("Trace [mG]");
	g_K1_Rm6p6->GetYaxis()->SetTitle("Trace [mG]");
	g_K3_R0->GetYaxis()->SetTitle("Trace [mG]");
	g_K1_R0_2022->GetYaxis()->SetTitle("Trace [mG]");

	g_K1_R0->SetLineWidth(2);
	g_K1_R3p2->SetLineWidth(2);
	g_K1_R6p6->SetLineWidth(2);
	g_K1_Rm6p6->SetLineWidth(2);
	g_K3_R0->SetLineWidth(2);
	g_K1_R0_2022->SetLineWidth(2);

	g_K1_R0->SetLineColor(1);
	g_K1_R3p2->SetLineColor(2);
	g_K1_R6p6->SetLineColor(3);
	g_K1_Rm6p6->SetLineColor(4);
	g_K3_R0->SetLineColor(5);
	g_K1_R0_2022->SetLineColor(6);


	int Npoints = 10080;
	double min_time = 1e9;
	double max_time = -1e9;

	for (int i=0; i<Npoints; i++){
		double time;
		double val;
		f0>>time;
		if (time<min_time) min_time = time;
		if (time>max_time) max_time = time;

		f1>>val;
		g_K1_R0->SetPoint(i,time,(val<-500?0:val));
		f2>>val;
		g_K1_R3p2->SetPoint(i,time,(val<-500?0:val));
		f3>>val;
		g_K1_R6p6->SetPoint(i,time,(val<-500?0:val));
		f4>>val;
		g_K1_Rm6p6->SetPoint(i,time,(val<-500?0:val));
		f5>>val;
		g_K3_R0->SetPoint(i,time,(val<-500?0:val));


		//f6>>time;
		f7>>val;
		g_K1_R0_2022->SetPoint(i,time,(val<-500?0:val));
	}

	double dt = (max_time-min_time)/(Npoints-1);

	TH1D* h1_K1_R0 = new TH1D("h1_K1_R0","K1 0.0 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);
	TH1D* h1_K1_R3p2 = new TH1D("h1_K1_R3p2","K1 3.2 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);
	TH1D* h1_K1_R6p6 = new TH1D("h1_K1_R6p6","K1 6.6 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);
	TH1D* h1_K1_Rm6p6 = new TH1D("h1_K1_Rm6p6","K1 -6.6 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);
	TH1D* h1_K3_R0 = new TH1D("h1_K3_R0","K3 0.0 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);
	TH1D* h1_K1_R0_2022 = new TH1D("h1_K1_R0_2022","K1 2022 0.0 mm",Npoints,min_time-0.5*dt,max_time+0.5*dt);

	for (int i=0; i<Npoints; i++){
		h1_K1_R0->Fill(g_K1_R0->GetPointX(i),g_K1_R0->GetPointY(i));
		h1_K1_R3p2->Fill(g_K1_R3p2->GetPointX(i),g_K1_R3p2->GetPointY(i));
		h1_K1_R6p6->Fill(g_K1_R6p6->GetPointX(i),g_K1_R6p6->GetPointY(i));
		h1_K1_Rm6p6->Fill(g_K1_Rm6p6->GetPointX(i),g_K1_Rm6p6->GetPointY(i));
		h1_K3_R0->Fill(g_K3_R0->GetPointX(i),g_K3_R0->GetPointY(i));
		h1_K1_R0_2022->Fill(g_K1_R0_2022->GetPointX(i),g_K1_R0_2022->GetPointY(i));
	}

	g_K1_R0->Write();
	g_K1_R3p2->Write();
	g_K1_R6p6->Write();
	g_K1_Rm6p6->Write();
	g_K3_R0->Write();
	g_K1_R0_2022->Write();

	h1_K1_R0->Write();
	h1_K1_R3p2->Write();
	h1_K1_R6p6->Write();
	h1_K1_Rm6p6->Write();
	h1_K3_R0->Write();
	h1_K1_R0_2022->Write();


	fout->Close();

}
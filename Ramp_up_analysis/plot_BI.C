void plot_BI(){


	TFile* fin = TFile::Open("output_Ramp_jan29_H25Q00_5173to0.root");
	TGraph* g_ramp = (TGraph*)fin->Get("Ramp_norm_current");

	vector<double> fit_pol2_centers = {
		156.8,
		791.5,
		1418.6,
		2037.6,
		2643.4,
		3280.0,
		3936.0,
		4697.4,
	};

	vector<double> fit_zero_centers = {
		461.7,
		1110.0,
		1740.5,
		2334.9,
		2947.8,
		3607.5,
		4281.8,
		5172.8
	};

	map<double,double> fit_points;

	new TCanvas();
	g_ramp->Draw("APL");

	TF1* f_pol2 = new TF1("f_pol2","[0]*(x-[1])*(x-[1])+[2]");
	for (int i=0; i<fit_pol2_centers.size(); i++){
		double center = fit_pol2_centers[i];
		double width = 100;
		f_pol2->SetParameters(0.1,center,g_ramp->Eval(center));
		TFitResultPtr fit = g_ramp->Fit(f_pol2,"S+","",center-width,center+width);
		fit_points[fit->Parameter(1)] = fit->ParError(1);
		f_pol2->DrawCopy("SAME");
	}

	TF1* f_pol1 = new TF1("f_pol1","[0]*(x-[1])");
	for (int i=0; i<fit_zero_centers.size(); i++){
		double center = fit_zero_centers[i];
		double width = 50;
		f_pol1->SetParameters(0.1,center);
		TFitResultPtr fit = g_ramp->Fit(f_pol1,"S+","",center-width,center+width);
		fit_points[fit->Parameter(1)] = fit->ParError(1);
		f_pol1->DrawCopy("SAME");
	}


	//Get sorted points
	vector<double> points;
	vector<double> points_errors;
	for (auto pair : fit_points){
		points.push_back(pair.first);
		points_errors.push_back(pair.second);
	}

	/*
	vector<double> points = {
		156.8,
		461.7,
		791.5,
		1110.0,
		1418.6,
		1740.5,
		2037.6,
		2334.9,
		2643.4,
		2947.8,
		3280.0,
		3607.5,
		3936.0,
		4281.8,
		4697.4,
		5172.8
	};

	vector<double> points_errors = {
		0.4,
		0.1,
		0.2,
		0.1,
		0.3,
		0.1,
		0.2,
		0.1,
		0.3,
		0.1,
		0.1,
		0.1,
		0.5,
		0.1,
		0.1,
		0.1
	};
	*/


	double Bnom = 1.45;
	double Ncycles = 15.5;

	TGraph* g_period = new TGraph();
	for (int i=1; i<points.size(); i++){
		double x = 0.5*(points[i-1]+points[i]);
		double y = points[i]-points[i-1];
		g_period->SetPoint(i-1,x,y);
	}

	new TCanvas();
	g_period->Draw("APL");

	TGraphErrors* g_BI = new TGraphErrors();

	for (int i=0; i<points.size(); i++){

		double B = Bnom - (points.size()-i-1)*(Bnom/Ncycles);
		g_BI->SetPoint(i,points[i],B);
		g_BI->SetPointError(i,points_errors[i],0);
	}
	g_BI->SetPoint(g_BI->GetN(),0,0);
	g_BI->Sort();


	new TCanvas("","",1000,1000);
	g_BI->SetMarkerStyle(20);
	g_BI->GetXaxis()->SetRangeUser(0,5300);
	g_BI->GetYaxis()->SetRangeUser(0,1.5);
	g_BI->GetXaxis()->SetTitle("Current [A]");
	g_BI->GetYaxis()->SetTitle("Bfield [T]");
	g_BI->Draw("AP");

	TLine* l1 = new TLine(0,0,5173,1.45);
	l1->SetLineWidth(2);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(kDashed);
	l1->Draw();

	TGraph* g_interpolation = new TGraph();
	for (double x=0; x<5173; x+=1){
		g_interpolation->AddPoint(x,g_BI->Eval(x,0,"S"));
	}
	g_interpolation->Draw("L");

	//TFile* fout = new TFile("BI.root","recreate");
	//g_BI->Write("g_BI");
	//g_interpolation->Write("g_BI_spline");
	//fout->Write();
	//fout->Close();

	TGraph* g_derivative = new TGraph();
	for (int i=1; i<g_BI->GetN(); i++){
		double x = g_BI->GetPointX(i);
		double prev_x = g_BI->GetPointX(i-1);
		double val = (g_BI->GetPointY(i)-g_BI->GetPointY(i-1))/(x-prev_x);
		g_derivative->AddPoint(0.5*(x+prev_x),1e4*val);
		cout<<0.5*(x+prev_x)<<" A -> "<<1e4*val<<" G/A\n";
	}
	new TCanvas();
	g_derivative->GetXaxis()->SetTitle("Current [A]");
	g_derivative->GetYaxis()->SetTitle("Field linearity [G/A]");
	g_derivative->SetMarkerStyle(20);
	g_derivative->Draw("APL");


}
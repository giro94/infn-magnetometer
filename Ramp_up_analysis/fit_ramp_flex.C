
void fit_ramp_flex(){

	vector<TString> filenames = {
		"output_Ramp_jan16_0to5175.root",
		"output_Ramp_jan16_5175to2000to5175.root",
		"output_Ramp_jan17_5175to0.root",
		"output_Ramp_jan18_0to5175.root",
		"output_Ramp_jan19_5173to4353.root",
		"output_Ramp_jan20_4353to3619.root",
		"output_Ramp_jan21_3619to3043.root",
		"output_Ramp_jan22_3043to5173.root",
		"output_Ramp_jan26_H25Q130_5173to2000to5173.root",
		"output_Ramp_jan29_H25Q00_5173to0.root"
	};

	double fit_range = 50; //Ampere
	double a_to_G = 1.45e4 / 5173.;
	double sensor_gain = 28.;

	int Nfiles = filenames.size();

	TFile** f = new TFile*[Nfiles];

	TGraphErrors** g_ramp = new TGraphErrors* [Nfiles];
	TGraphErrors** g_ramp_down = new TGraphErrors* [Nfiles];
	TGraphErrors** g_ramp_up = new TGraphErrors* [Nfiles];
	TGraph** g_rampdown_derivative = new TGraph* [Nfiles];
	TGraph** g_rampup_derivative = new TGraph* [Nfiles];
	vector<vector<double>> nodes_down;
	vector<vector<double>> nodes_up;
	nodes_down.resize(Nfiles);
	nodes_up.resize(Nfiles);
	vector<vector<double>> flex_down;
	vector<vector<double>> flex_up;
	flex_down.resize(Nfiles);
	flex_up.resize(Nfiles);


	TF1* f_pol1 = new TF1("f_pol1","[0]+[2]*(x-[1])*(x-[1])");

	for (int i=0; i<Nfiles; i++){

		f[i] = TFile::Open(filenames[i]);
		g_ramp[i] = (TGraphErrors*)f[i]->Get("Rampup");
		g_ramp_down[i] = new TGraphErrors();
		g_ramp_up[i] = new TGraphErrors();

		TString histTitle = filenames[i];
		histTitle.Remove(0,histTitle.Index("Ramp_")+5);
		g_ramp[i]->SetTitle(histTitle);
		g_ramp[i]->GetXaxis()->SetRangeUser(0,5500);
		g_ramp[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp[i]->SetLineWidth(2);
		g_ramp[i]->SetMarkerStyle(20);
		g_ramp[i]->SetMarkerColor(i%8+1);

		g_ramp_down[i]->SetTitle(histTitle);
		g_ramp_down[i]->GetXaxis()->SetRangeUser(0,5500);
		g_ramp_down[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp_down[i]->SetLineWidth(2);
		g_ramp_down[i]->SetMarkerStyle(20);
		g_ramp_down[i]->SetMarkerColor(i%8+1);
	
		g_ramp_up[i]->SetTitle(histTitle);
		g_ramp_up[i]->GetXaxis()->SetRangeUser(0,5500);
		g_ramp_up[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp_up[i]->SetLineWidth(2);
		g_ramp_up[i]->SetMarkerStyle(20);
		g_ramp_up[i]->SetMarkerColor(i%8+1);
	
	
		double minimum_current = 1e9;
		for (int j=0; j<g_ramp[i]->GetN(); j++){
			if (g_ramp[i]->GetPointX(j) < minimum_current){
				minimum_current = g_ramp[i]->GetPointX(j);
			}
		}

		bool down = true;
		for (int j=0; j<g_ramp[i]->GetN(); j++){
			if (g_ramp[i]->GetPointX(j) == minimum_current){
				down = false;
			}

			if (down){
				g_ramp_down[i]->SetPoint(g_ramp_down[i]->GetN(),g_ramp[i]->GetPointX(j),g_ramp[i]->GetPointY(j));
				g_ramp_down[i]->SetPointError(g_ramp_down[i]->GetN()-1,g_ramp[i]->GetErrorX(j),g_ramp[i]->GetErrorY(j));
			} else {
				g_ramp_up[i]->SetPoint(g_ramp_up[i]->GetN(),g_ramp[i]->GetPointX(j),g_ramp[i]->GetPointY(j));
				g_ramp_up[i]->SetPointError(g_ramp_up[i]->GetN()-1,g_ramp[i]->GetErrorX(j),g_ramp[i]->GetErrorY(j));
			}
		}


		g_rampdown_derivative[i] = new TGraph();
		g_rampup_derivative[i] = new TGraph();

		//Calculate derivatives
		for (int j=1; j<g_ramp_down[i]->GetN(); j++){
			double prev_x = g_ramp_down[i]->GetPointX(j-1);
			double prev_y = g_ramp_down[i]->GetPointY(j-1);
			double x = g_ramp_down[i]->GetPointX(j);
			double y = g_ramp_down[i]->GetPointY(j);
			double derivative = (y-prev_y)/(x-prev_x);
			//Suppress infinites
			if (abs(x-prev_x)<1) derivative = 0;
			if (abs(derivative) > 1) derivative = 0;

			g_rampdown_derivative[i]->SetPoint(j-1,0.5*(x+prev_x),derivative);
		}
		for (int j=1; j<g_ramp_up[i]->GetN(); j++){
			double prev_x = g_ramp_up[i]->GetPointX(j-1);
			double prev_y = g_ramp_up[i]->GetPointY(j-1);
			double x = g_ramp_up[i]->GetPointX(j);
			double y = g_ramp_up[i]->GetPointY(j);
			double derivative = (y-prev_y)/(x-prev_x);
			//Suppress infinites
			if (abs(x-prev_x)<1) derivative = 0;
			if (abs(derivative) > 1) derivative = 0;

			g_rampup_derivative[i]->SetPoint(j-1,0.5*(x+prev_x),derivative);
		}

		//Find original nodes
		for (int j=1; j<g_ramp_down[i]->GetN(); j++){
			double prev_x = g_ramp_down[i]->GetPointX(j-1);
			double prev_y = g_ramp_down[i]->GetPointY(j-1);
			double x = g_ramp_down[i]->GetPointX(j);
			double y = g_ramp_down[i]->GetPointY(j);
			if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
				double ratio = -prev_y/(y-prev_y);
				double x_node = ratio * (x-prev_x) + prev_x;
				nodes_down[i].push_back(x_node);
			}
		}
		sort(nodes_down[i].begin(),nodes_down[i].end());

		for (int j=1; j<g_ramp_up[i]->GetN(); j++){
			double prev_x = g_ramp_up[i]->GetPointX(j-1);
			double prev_y = g_ramp_up[i]->GetPointY(j-1);
			double x = g_ramp_up[i]->GetPointX(j);
			double y = g_ramp_up[i]->GetPointY(j);
			if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
				double ratio = -prev_y/(y-prev_y);
				double x_node = ratio * (x-prev_x) + prev_x;
				nodes_up[i].push_back(x_node);
			}
		}
		sort(nodes_up[i].begin(),nodes_up[i].end());

		//Find derivatives peaks/valleys
		for (int k=0; k<nodes_down[i].size(); k++){
			double d_min = 1e9;
			double d_max = -1e9;
			double x_min = -1;
			double x_max = -1;
			for (int j=1; j<g_rampdown_derivative[i]->GetN()-1; j++){
				double x = g_rampdown_derivative[i]->GetPointX(j);
				double y = g_rampdown_derivative[i]->GetPointY(j);
				if (x < nodes_down[i][k]-100 || x > nodes_down[i][k]+100) continue;

				if (y > d_max){
					d_max = y;
					x_max = x;
				}
				if (y < d_min){
					d_min = y;
					x_min = x;
				}
			}
			if (x_min < 0) flex_down[i].push_back(x_min);
			if (x_max > 0) flex_down[i].push_back(x_max);
		}
		sort(flex_down[i].begin(),flex_down[i].end());

		for (int k=0; k<nodes_up[i].size(); k++){
			double d_min = 1e9;
			double d_max = -1e9;
			double x_min = -1;
			double x_max = -1;
			for (int j=1; j<g_rampup_derivative[i]->GetN()-1; j++){
				double x = g_rampup_derivative[i]->GetPointX(j);
				double y = g_rampup_derivative[i]->GetPointY(j);
				if (x < nodes_up[i][k]-100 || x > nodes_up[i][k]+100) continue;

				if (y > d_max){
					d_max = y;
					x_max = x;
				}
				if (y < d_min){
					d_min = y;
					x_min = x;
				}
			}
			if (x_min < 0) flex_up[i].push_back(x_min);
			if (x_max > 0) flex_up[i].push_back(x_max);
		}
		sort(flex_up[i].begin(),flex_up[i].end());
	}




	gStyle->SetOptStat(0);
	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp[i]->Draw(i==0?"APL":"PL");
	}


	TGraphErrors** g_slope = new TGraphErrors* [Nfiles];
	TH1D* h1_slope = new TH1D("h1_slope","Slope distribution;Slope [mV/mG]",20,0.1,0.7);
	for (int i=0; i<Nfiles; i++){

		TCanvas* can = new TCanvas();
		can->Divide(3,1);
		can->cd(1);
		g_ramp[i]->Draw("APL");
		g_ramp_down[i]->Draw("PL");
		g_ramp_up[i]->Draw("PL");

		can->cd(2);
		g_rampdown_derivative[i]->GetXaxis()->SetLimits(0,5500);
		g_rampdown_derivative[i]->GetYaxis()->SetRangeUser(-0.1,0.1);
		g_rampup_derivative[i]->GetXaxis()->SetLimits(0,5500);
		g_rampup_derivative[i]->GetYaxis()->SetRangeUser(-0.1,0.1);
		g_rampdown_derivative[i]->Draw("APL");
		g_rampup_derivative[i]->Draw("PL");


		g_slope[i] = new TGraphErrors();
		for (int j=0; j<flex_down[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_down[i]->Fit("pol1","QS+","",flex_down[i][j]-fit_range,flex_down[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),flex_down[i][j],sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,sensor_gain*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill(sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
			}
		}

		for (int j=0; j<flex_up[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_up[i]->Fit("pol1","QS+","",flex_up[i][j]-fit_range,flex_up[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),flex_up[i][j],sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,sensor_gain*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill(sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
			}
		}
		g_slope[i]->Sort();
		
		can->cd(3);
		g_slope[i]->GetXaxis()->SetTitle("Current [A]");
		g_slope[i]->GetYaxis()->SetTitle("Slope [mV/mG]");
		g_slope[i]->GetXaxis()->SetLimits(0,5500);
		g_slope[i]->GetYaxis()->SetRangeUser(0.1,0.7);
		g_slope[i]->SetMarkerStyle(20);
		g_slope[i]->SetMarkerColor(i%8+1);
		g_slope[i]->Draw("APL");
	}


	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_slope[i]->Draw(i==0?"APLZ":"PLZ");
	}
	gPad->SetGridy();








	new TCanvas();
	h1_slope->Draw("HIST");

















}
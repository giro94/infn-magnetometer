
void fit_ramp_nodes(){

	vector<TString> filenames = {
		//"output_Ramp_jan16_0to5175.root",
		//"output_Ramp_jan16_5175to2000to5175.root",
		//"output_Ramp_jan17_5175to0.root",
		//"output_Ramp_jan18_0to5175.root",
		//"output_Ramp_jan19_5173to4353.root",
		//"output_Ramp_jan20_4353to3619.root",
		//"output_Ramp_jan21_3619to3043.root",
		//"output_Ramp_jan22_3043to5173.root",
		//"output_Ramp_jan26_H25Q130_5173to2000to5173.root",
		//"output_Ramp_jan29_H25Q00_5173to0.root"
		"output_FD_R0_ramp_oct9_H22p5.root",
		"output_FD_R1_ramp_oct8_H0.root"
	};

	bool use_normalized = true;
	bool calibrate = true;

	TFile* f_calibration = TFile::Open("BI.root");
	TGraphErrors* g_BI = (TGraphErrors*)f_calibration->Get("g_BI");

	double fit_range = calibrate?0.014:50; //Ampere or teslas
	double a_to_G = calibrate?1e4:(1.45e4 / 5173.);
	double sensor_gain = 28.13;

	int Nfiles = filenames.size();

	TFile** f = new TFile*[Nfiles];

	TGraphErrors** g_ramp = new TGraphErrors* [Nfiles];
	TGraphErrors** g_ramp_down = new TGraphErrors* [Nfiles];
	TGraphErrors** g_ramp_up = new TGraphErrors* [Nfiles];
	vector<vector<double>> nodes_down;
	vector<vector<double>> nodes_up;
	nodes_down.resize(Nfiles);
	nodes_up.resize(Nfiles);



	for (int i=0; i<Nfiles; i++){

		cout<<filenames[i]<<"\n";
		f[i] = TFile::Open(filenames[i]);
		g_ramp[i] = (TGraphErrors*)f[i]->Get(use_normalized?"Ramp_norm":"Ramp");
		if (g_ramp[i] == nullptr) g_ramp[i] = (TGraphErrors*)f[i]->Get("Rampup");
		g_ramp_down[i] = new TGraphErrors();
		g_ramp_up[i] = new TGraphErrors();

		if (calibrate){
			for (int j=0; j<g_ramp[i]->GetN(); j++){
				g_ramp[i]->SetPointX(j,g_BI->Eval(g_ramp[i]->GetPointX(j),0,"S"));
				g_ramp[i]->GetXaxis()->SetTitle("Bfield [T]");
			}
		}

		//hack to normalize H25Q00 to H25Q130 (blumlein ratio = 0.69)
		//if (filenames[i]=="output_Ramp_jan29_H25Q00_5173to0.root"){
		//	for (int j=0; j<g_ramp[i]->GetN(); j++){
		//		g_ramp[i]->SetPoint(j,g_ramp[i]->GetPointX(j),g_ramp[i]->GetPointY(j)/0.69);
		//		g_ramp[i]->SetPointError(j,g_ramp[i]->GetErrorX(j),g_ramp[i]->GetErrorY(j)/0.69);
		//	}
		//}

		TString histTitle = filenames[i];
		histTitle.Remove(0,histTitle.Index("Ramp_")+5);
		g_ramp[i]->SetTitle(histTitle);
		g_ramp[i]->GetXaxis()->SetRangeUser(0,calibrate?1.5:5500);
		g_ramp[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp[i]->SetLineWidth(2);
		g_ramp[i]->SetMarkerStyle(20);
		g_ramp[i]->SetMarkerColor(i%8+1);

		g_ramp_down[i]->SetTitle(histTitle);
		g_ramp_down[i]->GetXaxis()->SetRangeUser(0,calibrate?1.5:5500);
		g_ramp_down[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp_down[i]->SetLineWidth(2);
		g_ramp_down[i]->SetMarkerStyle(20);
		g_ramp_down[i]->SetMarkerColor(i%8+1);
	
		g_ramp_up[i]->SetTitle(histTitle);
		g_ramp_up[i]->GetXaxis()->SetRangeUser(0,calibrate?1.5:5500);
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

		for (int j=1; j<g_ramp_down[i]->GetN(); j++){
			double prev_x = g_ramp_down[i]->GetPointX(j-1);
			double prev_y = g_ramp_down[i]->GetPointY(j-1);
			double x = g_ramp_down[i]->GetPointX(j);
			double y = g_ramp_down[i]->GetPointY(j);
			if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
				double ratio = -prev_y/(y-prev_y);
				double x_node = ratio * (x-prev_x) + prev_x;
				if (abs(x-prev_x)>(calibrate?3e-5:0.1)){
					nodes_down[i].push_back(x_node);
				}
			}
		}
		//Artificially put last node
		nodes_down[i].push_back(calibrate?1.438:5130);
		sort(nodes_down[i].begin(),nodes_down[i].end());

		for (int j=1; j<g_ramp_up[i]->GetN(); j++){
			double prev_x = g_ramp_up[i]->GetPointX(j-1);
			double prev_y = g_ramp_up[i]->GetPointY(j-1);
			double x = g_ramp_up[i]->GetPointX(j);
			double y = g_ramp_up[i]->GetPointY(j);
			if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
				double ratio = -prev_y/(y-prev_y);
				double x_node = ratio * (x-prev_x) + prev_x;
				if (abs(x-prev_x)>(calibrate?3e-5:0.1)){
					nodes_up[i].push_back(x_node);
				}
			}
		}
		sort(nodes_up[i].begin(),nodes_up[i].end());

		f[i]->Close();
	}


	gStyle->SetOptStat(0);


	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp[i]->Draw(i==0?"APL":"PL");
	}

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp_down[i]->Draw(i==0?"APL":"PL");
	}
	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp_up[i]->Draw(i==0?"APL":"PL");
	}

	TGraphErrors** g_slope = new TGraphErrors* [Nfiles];
	TH1D* h1_slope = new TH1D("h1_slope","Slope distribution;Slope [mV/mG]",20,0.1,0.7);
	for (int i=0; i<Nfiles; i++){
		g_slope[i] = new TGraphErrors();

		new TCanvas();
		g_ramp[i]->Draw("APL");
		g_ramp_down[i]->Draw("PL");
		g_ramp_up[i]->Draw("PL");

		for (int j=0; j<nodes_down[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_down[i]->Fit("pol1","QS+","",nodes_down[i][j]-fit_range,nodes_down[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),nodes_down[i][j],sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,sensor_gain*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill(sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
			}
		}

		for (int j=0; j<nodes_up[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_up[i]->Fit("pol1","QS+","",nodes_up[i][j]-fit_range,nodes_up[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),nodes_up[i][j],sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,sensor_gain*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill(sensor_gain*abs(fit_res->Parameter(1)/a_to_G));
			}
		}
	}


	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_slope[i]->Sort();
		
		g_slope[i]->SetName(Form("g_slope_%d",i));
		g_slope[i]->SetTitle(Form("g_slope_%d",i));
		g_slope[i]->GetXaxis()->SetTitle(calibrate?"Bfield [T]":"Current [A]");
		g_slope[i]->GetYaxis()->SetTitle("Slope [mV/mG]");
		g_slope[i]->GetXaxis()->SetLimits(0,calibrate?1.5:5500);
		g_slope[i]->GetYaxis()->SetRangeUser(0.1,0.7);
		g_slope[i]->SetMarkerStyle(20);
		g_slope[i]->SetMarkerColor(i%8+1);
		g_slope[i]->Draw(i==0?"APLZ":"PLZ");

		TFitResultPtr res = g_slope[i]->Fit("pol0","SQ+","",0,calibrate?1.4:4000);
		cout<<"Slope fit: "<<res->Parameter(0)<<" +- "<<res->ParError(0)<<"\n";
		double y=0;
		double y2=0;
		double nfit=0;
		for (int j=0; j<g_slope[i]->GetN(); j++){
			if (g_slope[i]->GetPointX(j) < (calibrate?1.5:4000)){
				y += g_slope[i]->GetPointY(j);
				y2 += g_slope[i]->GetPointY(j)*g_slope[i]->GetPointY(j);
				nfit += 1;
			}
		}
		y /= nfit;
		y2 /= nfit;
		double rms = sqrt(y2 - y*y);
		cout<<"RMS: "<<rms<<", mean error: "<<rms/sqrt(nfit)<<"\n";
	}
	gPad->SetGridy();


	new TCanvas();
	h1_slope->Draw("HIST");




















}
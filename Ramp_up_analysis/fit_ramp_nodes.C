
void fit_ramp_nodes(){

	vector<TString> filenames = {
	//	"output_Rampup_R0_H20_oct11.root",
	//	"output_Rampup_R0_H20_oct13.root",
		//"output_Rampup_R0_H25_oct20_B50-88.root",
		//"output_Rampup_R0_H25_oct20.root",
		//"output_Rampup_R0_H25_oct21_B88-100.root",
		//"output_Rampup_R0_H25_oct24_B0-100.root",
	//	"output_Rampup_R0_H30_oct16.root",
	//	"output_Rampup_R1_H5_oct5.root",

		//"output_Ramp_jan16_0to5175.root",
		//"output_Ramp_jan16_5175to2000to5175.root",
		//"output_Ramp_jan17_5175to0.root",
		//"output_Ramp_jan18_0to5175.root",
		//"output_Ramp_jan19_5173to4353.root",
		//"output_Ramp_jan20_4353to3619.root",
		//"output_Ramp_jan21_3619to3043.root",
		//"output_Ramp_jan22_3043to5173.root",
		//"output_Ramp_jan26_H25Q130_5173to2000to5173.root",
		"output_Ramp_jan29_H25Q00_5173to0.root",
		
		//"output_FD_R0_ramp_oct9_H22p5.root",
		//"output_FD_R1_ramp_oct8_H0.root"
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
		g_ramp[i] = (TGraphErrors*)f[i]->Get(use_normalized?"Ramp_norm_current":"Ramp_current");
		if (g_ramp[i] == nullptr) g_ramp[i] = (TGraphErrors*)f[i]->Get("Rampup_current");
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
		histTitle.Remove(0,histTitle.Index("Ramp")+5);
		g_ramp[i]->SetTitle(histTitle);
		g_ramp[i]->GetXaxis()->SetRangeUser(0,calibrate?1.5:5500);
		g_ramp[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp[i]->SetLineWidth(2);
		g_ramp[i]->SetMarkerStyle(20);
		g_ramp[i]->SetMarkerColor(i%8+1);

		g_ramp_down[i]->SetTitle(histTitle+" (down)");
		g_ramp_down[i]->GetXaxis()->SetRangeUser(0,calibrate?1.5:5500);
		g_ramp_down[i]->GetYaxis()->SetRangeUser(-15,15);
		g_ramp_down[i]->SetLineWidth(2);
		g_ramp_down[i]->SetMarkerStyle(20);
		g_ramp_down[i]->SetMarkerColor(i%8+1);
	
		g_ramp_up[i]->SetTitle(histTitle+" (up)");
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


	double sine_th_min = 1;
	double sine_th_max = 12;
	double sine_th_step = 0.1;
	int Nsteps = (sine_th_max-sine_th_min)/sine_th_step;

	TGraphErrors** g_slope = new TGraphErrors* [Nfiles];
	TGraphErrors*** g_ramp_sine = new TGraphErrors** [Nfiles];
	TGraphErrors** g_slope_sine = new TGraphErrors* [Nfiles];
	TGraphErrors** g_chi2_sine = new TGraphErrors* [Nfiles];
	TH1D* h1_slope = new TH1D("h1_slope","Slope distribution;Slope [mV/mG]",20,0.1,0.7);
	TF1* f_sine = new TF1("f_sine","[0]*sin([1]*x+[2])",0,1.45);
	for (int i=0; i<Nfiles; i++){
		g_slope[i] = new TGraphErrors();

		bool fast_diodes = (filenames[i].Contains("FD"));
		
		new TCanvas();
		g_ramp[i]->Draw("APL");
		g_ramp_down[i]->Draw("PL");
		g_ramp_up[i]->Draw("PL");

		// Fit nodes rampdown
		for (int j=0; j<nodes_down[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_down[i]->Fit("pol1","QS+","",nodes_down[i][j]-fit_range,nodes_down[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),nodes_down[i][j],(fast_diodes?1:sensor_gain)*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,(fast_diodes?1:sensor_gain)*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill((fast_diodes?1:sensor_gain)*abs(fit_res->Parameter(1)/a_to_G));
			}
		}

		// Fit nodes rampup
		for (int j=0; j<nodes_up[i].size(); j++){
			TFitResultPtr fit_res = g_ramp_up[i]->Fit("pol1","QS+","",nodes_up[i][j]-fit_range,nodes_up[i][j]+fit_range);
			if (fit_res >= 0){
				g_slope[i]->SetPoint(g_slope[i]->GetN(),nodes_up[i][j],(fast_diodes?1:sensor_gain)*abs(fit_res->Parameter(1)/a_to_G));
				g_slope[i]->SetPointError(g_slope[i]->GetN()-1,0,(fast_diodes?1:sensor_gain)*fit_res->ParError(1)/a_to_G);
				h1_slope->Fill((fast_diodes?1:sensor_gain)*abs(fit_res->Parameter(1)/a_to_G));
			}
		}

		// Fit sine rampdown
		g_ramp_sine[i] = new TGraphErrors* [Nsteps];
		g_slope_sine[i] = new TGraphErrors();
		g_chi2_sine[i] = new TGraphErrors();
		for (int j=0; j<Nsteps; j++){
			g_ramp_sine[i][j] = new TGraphErrors();
			double sine_th = sine_th_min+j*sine_th_step;
			for (int k=0; k<g_ramp_up[i]->GetN(); k++){
				if (abs(g_ramp_up[i]->GetPointY(k)) < sine_th){
					g_ramp_sine[i][j]->SetPoint(g_ramp_sine[i][j]->GetN(),g_ramp_up[i]->GetPointX(k),g_ramp_up[i]->GetPointY(k));
				}
			}
			f_sine->SetParameters(10,8*M_PI/1.45,0);
			TFitResultPtr fit_res = g_ramp_sine[i][j]->Fit("f_sine","QS+","",0,1.45);
			if (fit_res >= 0){
				double slope_fit = (fast_diodes?1:sensor_gain)*abs(fit_res->Parameter(0)*fit_res->Parameter(1)/a_to_G);
				double slope_fit_err = slope_fit*(abs(fit_res->ParError(0)/fit_res->Parameter(0)) + abs(fit_res->ParError(1)/fit_res->Parameter(1)));
				g_slope_sine[i]->SetPoint(g_slope_sine[i]->GetN(),sine_th,slope_fit);
				g_slope_sine[i]->SetPointError(g_slope_sine[i]->GetN()-1,0,slope_fit_err);
				g_chi2_sine[i]->SetPoint(g_chi2_sine[i]->GetN(),sine_th,fit_res->Chi2()/fit_res->Ndf());
			}
		}
	}


	new TCanvas();
	TLegend* leg_fit = new TLegend(0.5,0.7,0.9,0.9);
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

		TFitResultPtr res = g_slope[i]->Fit("pol0","SQ+","",calibrate?0.1:285,calibrate?1.4:4000);
		//TFitResultPtr res = g_slope[i]->Fit("pol0","SQ+","",1,1.4);
		cout<<"Slope fit: "<<res->Parameter(0)<<" +- "<<res->ParError(0)<<"\n";
		double y=0;
		double y2=0;
		double nfit=0;
		for (int j=0; j<g_slope[i]->GetN(); j++){
			if (g_slope[i]->GetPointX(j) < (calibrate?1.5:4000) && g_slope[i]->GetPointX(j) > (calibrate?0.1:285)){
				y += g_slope[i]->GetPointY(j);
				y2 += g_slope[i]->GetPointY(j)*g_slope[i]->GetPointY(j);
				nfit += 1;
			}
		}
		y /= nfit;
		y2 /= nfit;
		double rms = sqrt(y2 - y*y);
		cout<<"RMS: "<<rms<<", mean error: "<<rms/sqrt(nfit)<<"\n";

		TString histTitle = filenames[i];
		histTitle.Remove(0,histTitle.Index("Ramp")+5);
		leg_fit->AddEntry(g_slope[i],Form("%s | %.3f +- %.3f mV/mG",histTitle.Data(),res->Parameter(0),rms/sqrt(nfit)),"PL");
	}
	leg_fit->Draw();
	gPad->SetGridy();


	new TCanvas();
	h1_slope->Draw("HIST");

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_slope_sine[i]->SetName(Form("g_slope_sine_%d",i));
		g_slope_sine[i]->SetTitle(Form("g_slope_sine_%d",i));
		g_slope_sine[i]->GetXaxis()->SetTitle("Cutoff [V]");
		g_slope_sine[i]->GetYaxis()->SetTitle("Slope [mV/mG]");
		g_slope_sine[i]->Draw(i==0?"APLZ":"PLZ");
	}

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_chi2_sine[i]->SetName(Form("g_chi2_sine_%d",i));
		g_chi2_sine[i]->SetTitle(Form("g_chi2_sine_%d",i));
		g_chi2_sine[i]->GetXaxis()->SetTitle("Cutoff [V]");
		g_chi2_sine[i]->GetYaxis()->SetTitle("Fit Chi2/NDF");
		g_chi2_sine[i]->Draw(i==0?"APLZ":"PLZ");
	}

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp_sine[i][0]->SetName(Form("g_ramp_sine_%d_%d",i,Nsteps/2));
		g_ramp_sine[i][0]->SetTitle(Form("g_ramp_sine_%d (Cutoff %.1f)",i,sine_th_min+0*sine_th_step));
		g_ramp_sine[i][0]->GetXaxis()->SetTitle(calibrate?"Bfield [T]":"Current [A]");
		g_ramp_sine[i][0]->GetYaxis()->SetTitle("B-A [V] (12 V)");
		g_ramp_sine[i][0]->SetMarkerStyle(20);
		g_ramp_sine[i][0]->SetMarkerColor(i%8+1);
		g_ramp_sine[i][0]->GetYaxis()->SetRangeUser(-12,12);
		g_ramp_sine[i][0]->Draw(i==0?"APLZ":"PLZ");
	}

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp_sine[i][Nsteps/2]->SetName(Form("g_ramp_sine_%d_%d",i,Nsteps/2));
		g_ramp_sine[i][Nsteps/2]->SetTitle(Form("g_ramp_sine_%d (Cutoff %.1f)",i,sine_th_min+Nsteps/2*sine_th_step));
		g_ramp_sine[i][Nsteps/2]->GetXaxis()->SetTitle(calibrate?"Bfield [T]":"Current [A]");
		g_ramp_sine[i][Nsteps/2]->GetYaxis()->SetTitle("B-A [V] (12 V)");
		g_ramp_sine[i][Nsteps/2]->SetMarkerStyle(20);
		g_ramp_sine[i][Nsteps/2]->SetMarkerColor(i%8+1);
		g_ramp_sine[i][Nsteps/2]->GetYaxis()->SetRangeUser(-12,12);
		g_ramp_sine[i][Nsteps/2]->Draw(i==0?"APLZ":"PLZ");
	}

	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		g_ramp_sine[i][Nsteps-1]->SetName(Form("g_ramp_sine_%d_%d",i,Nsteps-1));
		g_ramp_sine[i][Nsteps-1]->SetTitle(Form("g_ramp_sine_%d (Cutoff %.1f)",i,sine_th_min+Nsteps-1*sine_th_step));
		g_ramp_sine[i][Nsteps-1]->GetXaxis()->SetTitle(calibrate?"Bfield [T]":"Current [A]");
		g_ramp_sine[i][Nsteps-1]->GetYaxis()->SetTitle("B-A [V] (12 V)");
		g_ramp_sine[i][Nsteps-1]->SetMarkerStyle(20);
		g_ramp_sine[i][Nsteps-1]->SetMarkerColor(i%8+1);
		g_ramp_sine[i][Nsteps-1]->GetYaxis()->SetRangeUser(-12,12);
		g_ramp_sine[i][Nsteps-1]->Draw(i==0?"APLZ":"PLZ");
	}

















}
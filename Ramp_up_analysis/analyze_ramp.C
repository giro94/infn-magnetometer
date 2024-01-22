

void analyze_ramp(TString filename){


	TFile* f1 = TFile::Open(filename);
	TGraph* g_A = (TGraph*)f1->Get("A");
	TGraph* g_B = (TGraph*)f1->Get("B");
	TGraph* g_ramp = (TGraph*)f1->Get("Rampup");
	TGraph* g_rampdown = new TGraph();
	TGraph* g_rampup = new TGraph();

	double minimum_current = 1e9;
	for (int i=0; i<g_ramp->GetN(); i++){
		if (g_ramp->GetPointX(i) < minimum_current){
			minimum_current = g_ramp->GetPointX(i);
		}
	}

	bool down = true;
	for (int i=0; i<g_ramp->GetN(); i++){
		if (g_ramp->GetPointX(i) == minimum_current){
			down = false;
		}

		if (down){
			g_rampdown->SetPoint(g_rampdown->GetN(),g_ramp->GetPointX(i),g_ramp->GetPointY(i));
		} else {
			g_rampup->SetPoint(g_rampup->GetN(),g_ramp->GetPointX(i),g_ramp->GetPointY(i));
		}
	}

	new TCanvas();
	g_rampdown->SetTitle("Rampdown;Current [A];B-A [V]");
	g_rampup->SetTitle("Rampup;Current [A];B-A [V]");
	g_rampdown->SetMarkerStyle(20);
	g_rampdown->SetMarkerColor(kBlue);
	g_rampup->SetMarkerStyle(20);
	g_rampup->SetMarkerColor(kRed);
	g_rampdown->Draw("APL");
	g_rampup->Draw("PL");

	vector<double> nodes_down;
	for (int i=1; i<g_rampdown->GetN(); i++){
		double prev_x = g_rampdown->GetPointX(i-1);
		double prev_y = g_rampdown->GetPointY(i-1);
		double x = g_rampdown->GetPointX(i);
		double y = g_rampdown->GetPointY(i);
		if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
			double ratio = -prev_y/(y-prev_y);
			double x_node = ratio * (x-prev_x) + prev_x;
			nodes_down.push_back(x_node);
		}
	}
	sort(nodes_down.begin(),nodes_down.end());

	vector<double> nodes_up;
	for (int i=1; i<g_rampup->GetN(); i++){
		double prev_x = g_rampup->GetPointX(i-1);
		double prev_y = g_rampup->GetPointY(i-1);
		double x = g_rampup->GetPointX(i);
		double y = g_rampup->GetPointY(i);
		if ((prev_y < 0 && y > 0) || (prev_y > 0 && y < 0)){
			double ratio = -prev_y/(y-prev_y);
			double x_node = ratio * (x-prev_x) + prev_x;
			nodes_up.push_back(x_node);
		}
	}
	sort(nodes_up.begin(),nodes_up.end());


	TGraph* g_ramp_derivative = new TGraph();
	TGraph* g_rampdown_derivative = new TGraph();
	TGraph* g_rampup_derivative = new TGraph();

	for (int i=1; i<g_ramp->GetN(); i++){
		double prev_x = g_ramp->GetPointX(i-1);
		double prev_y = g_ramp->GetPointY(i-1);
		double x = g_ramp->GetPointX(i);
		double y = g_ramp->GetPointY(i);
		double derivative = (y-prev_y)/(x-prev_x);
		//Suppress infinites
		if (abs(x-prev_x)<1) derivative = 0;

		g_ramp_derivative->SetPoint(i-1,0.5*(x+prev_x),50*derivative);
	}
	for (int i=1; i<g_rampdown->GetN(); i++){
		double prev_x = g_rampdown->GetPointX(i-1);
		double prev_y = g_rampdown->GetPointY(i-1);
		double x = g_rampdown->GetPointX(i);
		double y = g_rampdown->GetPointY(i);
		double derivative = (y-prev_y)/(x-prev_x);
		//Suppress infinites
		if (abs(x-prev_x)<1) derivative = 0;

		g_rampdown_derivative->SetPoint(i-1,0.5*(x+prev_x),50*derivative);
	}
	for (int i=1; i<g_rampup->GetN(); i++){
		double prev_x = g_rampup->GetPointX(i-1);
		double prev_y = g_rampup->GetPointY(i-1);
		double x = g_rampup->GetPointX(i);
		double y = g_rampup->GetPointY(i);
		double derivative = (y-prev_y)/(x-prev_x);
		//Suppress infinites
		if (abs(x-prev_x)<1) derivative = 0;

		g_rampup_derivative->SetPoint(i-1,0.5*(x+prev_x),50*derivative);
	}

	new TCanvas();
	g_ramp->Draw("APL");
	g_ramp_derivative->SetTitle("Derivative");
	g_ramp_derivative->GetXaxis()->SetTitle("Current [A]");
	g_ramp_derivative->GetYaxis()->SetTitle("Derivative [V/A]");
	g_ramp_derivative->SetMarkerStyle(20);
	g_ramp_derivative->SetMarkerColor(kRed);
	g_ramp_derivative->Draw("L");


	new TCanvas();
	g_rampdown->Draw("APL");
	g_rampdown_derivative->SetTitle("Derivative");
	g_rampdown_derivative->GetXaxis()->SetTitle("Current [A]");
	g_rampdown_derivative->GetYaxis()->SetTitle("Derivative [V/A]");
	g_rampdown_derivative->SetMarkerStyle(20);
	g_rampdown_derivative->SetMarkerColor(kBlue);
	g_rampdown_derivative->SetLineWidth(2);
	g_rampdown_derivative->Draw("L");

	TF1* f_pol1 = new TF1("f_pol1","[0]+[2]*(x-[1])*(x-[1])");

	vector<double> guesses_down = {2400,3000,3600,4300};
	vector<double> fits_down;
	f_pol1->SetParameters(2,guesses_down[0],-0.01);
	g_rampdown_derivative->Fit("f_pol1","Q+","",guesses_down[0]-150,guesses_down[0]+200);
	fits_down.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_down[1],-0.01);
	g_rampdown_derivative->Fit("f_pol1","Q+","",guesses_down[1]-150,guesses_down[1]+200);
	fits_down.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_down[2],-0.01);
	g_rampdown_derivative->Fit("f_pol1","Q+","",guesses_down[2]-150,guesses_down[2]+200);
	fits_down.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_down[3],-0.01);
	g_rampdown_derivative->Fit("f_pol1","Q+","",guesses_down[3]-150,guesses_down[3]+200);
	fits_down.push_back(f_pol1->GetParameter(1));

	TGraph* g_down_fits = new TGraph();
	for (int i=0; i<4; i++){
		double y_fit = g_rampdown->Eval(fits_down[i]);
		g_down_fits->SetPoint(i,fits_down[i],y_fit);
	}
	g_down_fits->SetMarkerStyle(89);
	g_down_fits->SetMarkerSize(2);
	g_down_fits->SetMarkerColor(kGreen);
	g_down_fits->Draw("P");
	TLine* line_zero = new TLine(minimum_current,0,5300,0);
	line_zero->Draw("L");
	gPad->SetGridy();

	TGraph* g_down_nodes = new TGraph();
	for (int i=0; i<4; i++){
		double y_node = g_rampdown->Eval(nodes_down[i]);
		g_down_nodes->SetPoint(i,nodes_down[i],y_node);
	}
	g_down_nodes->SetMarkerStyle(89);
	g_down_nodes->SetMarkerSize(2);
	g_down_nodes->SetMarkerColor(kGreen);
	g_down_nodes->Draw("P");


	new TCanvas();
	g_rampup->Draw("APL");
	g_rampup_derivative->SetTitle("Derivative");
	g_rampup_derivative->GetXaxis()->SetTitle("Current [A]");
	g_rampup_derivative->GetYaxis()->SetTitle("Derivative [V/A]");
	g_rampup_derivative->SetMarkerStyle(20);
	g_rampup_derivative->SetMarkerColor(kRed);
	g_rampup_derivative->SetLineWidth(2);
	g_rampup_derivative->Draw("L");

	vector<double> guesses_up = {2450,3100,3700,4400};
	vector<double> fits_up;
	f_pol1->SetParameters(2,guesses_up[0],-0.01);
	g_rampup_derivative->Fit("f_pol1","Q+","",guesses_up[0]-150,guesses_up[0]+200);
	fits_up.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_up[1],-0.01);
	g_rampup_derivative->Fit("f_pol1","Q+","",guesses_up[1]-50,guesses_up[1]+100);
	fits_up.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_up[2],-0.01);
	g_rampup_derivative->Fit("f_pol1","Q+","",guesses_up[2]-150,guesses_up[2]+200);
	fits_up.push_back(f_pol1->GetParameter(1));
	f_pol1->SetParameters(2,guesses_up[3],-0.01);
	g_rampup_derivative->Fit("f_pol1","Q+","",guesses_up[3]-200,guesses_up[3]+50);
	fits_up.push_back(f_pol1->GetParameter(1));

	TGraph* g_up_fits = new TGraph();
	for (int i=0; i<4; i++){
		double y_fit = g_rampup->Eval(fits_up[i]);
		g_up_fits->SetPoint(i,fits_up[i],y_fit);
	}
	g_up_fits->SetMarkerStyle(89);
	g_up_fits->SetMarkerSize(2);
	g_up_fits->SetMarkerColor(kGreen);
	g_up_fits->Draw("P");
	line_zero->Draw("L");
	gPad->SetGridy();

	TGraph* g_up_nodes = new TGraph();
	for (int i=0; i<4; i++){
		double y_node = g_rampup->Eval(nodes_up[i]);
		g_up_nodes->SetPoint(i,nodes_up[i],y_node);
	}
	g_up_nodes->SetMarkerStyle(89);
	g_up_nodes->SetMarkerSize(2);
	g_up_nodes->SetMarkerColor(kGreen);
	g_up_nodes->Draw("P");

	cout<<"Setpoints ramp down zero-crossing:\n";
	for (int i=0; i<4; i++){
		cout<<round(nodes_down[i])<<" ";
	}
	cout<<"\n";

	cout<<"Setpoints ramp down derivative:\n";
	for (int i=0; i<4; i++){
		cout<<round(fits_down[i])<<" ";
	}
	cout<<"\n";

	cout<<"Setpoints ramp up zero-crossing:\n";
	for (int i=0; i<4; i++){
		cout<<round(nodes_up[i])<<" ";
	}
	cout<<"\n";

	cout<<"Setpoints ramp up derivative:\n";
	for (int i=0; i<4; i++){
		cout<<round(fits_up[i])<<" ";
	}
	cout<<"\n";
}
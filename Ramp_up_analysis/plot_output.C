

void plot_output(TString filename){

	TFile* f_calibration = TFile::Open("BI.root");
	TGraphErrors* g_BI = (TGraphErrors*)f_calibration->Get("g_BI");

	TFile* f_in = TFile::Open(filename);

	TGraph* g_ramp = (TGraph*)f_in->Get("Ramp");
	TGraph* g_ramp_norm = (TGraph*)f_in->Get("Ramp_norm");

	TGraph* g_A = (TGraph*)f_in->Get("A");
	TGraph* g_B = (TGraph*)f_in->Get("B");
	TGraph* g_AB = (TGraph*)f_in->Get("AB");
	g_AB->SetTitle("A+B [V]");


	TGraph* g_ramp_calib = (TGraph*)g_ramp->Clone("g_ramp_calib");
	TGraph* g_ramp_norm_calib = (TGraph*)g_ramp_norm->Clone("g_ramp_norm_calib");

	for (int i=0; i<g_ramp->GetN(); i++){
		g_ramp_calib->SetPoint(i,g_BI->Eval(g_ramp->GetPointX(i),0,"S"),g_ramp->GetPointY(i));
		g_ramp_calib->GetXaxis()->SetTitle("Bfield [T]");

		g_ramp_norm_calib->SetPoint(i,g_BI->Eval(g_ramp_norm->GetPointX(i),0,"S"),g_ramp_norm->GetPointY(i));
		g_ramp_norm_calib->GetXaxis()->SetTitle("Bfield [T]");
	}



	new TCanvas("","",1000,800);
	g_A->Draw("APL");
	g_B->Draw("PL");
	g_AB->Draw("PL");
	gPad->SetGridy();
	gPad->BuildLegend();

	new TCanvas("","",1000,800);
	g_ramp_norm->SetMarkerColor(8);
	g_ramp->Draw("APL");
	g_ramp_norm->Draw("PL");
	gPad->SetGridy();
	gPad->BuildLegend();


	new TCanvas("","",1000,800);
	g_ramp_norm_calib->SetMarkerColor(8);
	g_ramp_calib->Draw("APL");
	g_ramp_norm_calib->Draw("PL");
	gPad->SetGridy();
	gPad->BuildLegend();







}
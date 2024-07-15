{


	TGraphErrors* g_blumlein = new TGraphErrors();
	TGraphErrors* g_transient = new TGraphErrors();
	g_blumlein->SetPoint(0,0,13.66);
	g_blumlein->SetPointError(0,0,0.03);
	g_blumlein->SetPoint(1,17.5,19.911);
	g_blumlein->SetPointError(1,0,0.008);
	g_blumlein->SetPoint(2,0,16.6);
	g_blumlein->SetPointError(2,0,0.39);


	g_transient->SetPoint(0,0,3.96);
	g_transient->SetPointError(0,0,0.06);
	g_transient->SetPoint(1,17.5,8.89);
	g_transient->SetPointError(1,0,0.03);


	g_blumlein->SetMarkerStyle(20);
	g_blumlein->SetMarkerColor(kBlack);

	g_transient->SetMarkerStyle(22);
	g_transient->SetMarkerColor(kRed);


	new TCanvas();
	g_blumlein->GetXaxis()->SetLimits(-20,20);
	g_blumlein->Draw("AP");
	g_transient->Draw("P");
}
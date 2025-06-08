{
	TFile* fin = TFile::Open("paper_radial_plot.root");
	TGraphErrors* g_all = (TGraphErrors*)fin->Get("g_all");
	TGraphErrors* g_INFN = (TGraphErrors*)fin->Get("g_INFN");
	TGraphErrors* g_UMass = (TGraphErrors*)fin->Get("g_UMass");
	TH1D* h1_model = (TH1D*)fin->Get("model");
	TGraphErrors* g_normpoint = (TGraphErrors*)fin->Get("normalization_point");

	TCanvas* can = new TCanvas("","",800,800);
	g_all->Draw("APZ");
	h1_umass_y0_INFN->Draw("HIST L SAME");
	g_INFN->Draw("PZ");
	g_UMass->Draw("PZ");
	gPad->SetGridx();
	gPad->SetGridy();
	TLegend* leg2 = new TLegend(0.35,0.6,0.65,0.8);
	leg2->AddEntry(g_INFN,"INFN","PL");
	leg2->AddEntry(g_UMass,"UMass","PL");
	leg2->AddEntry(h1_umass_y0_INFN,"UMass radial model","L");
	leg2->Draw();
	g_normpoint->Draw("P");

}
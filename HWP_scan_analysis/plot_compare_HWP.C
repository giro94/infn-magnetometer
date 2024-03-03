void plot_compare_HWP(){

	vector<TString> filenames = {
		"HWPscan_jan21_B3043_output.root",
		"HWPscan_jan20_B3619_output.root",
		"HWPscan_jan19_B4353_output.root",
		"HWPscan_jan22_B5173_output.root",
		//"HWPscan_jan23_B5173_output.root"
	};

	//vector<TString> filenames = {
	//	"HWPscan_jan16_B5175_output.root",
	//	"HWPscan_jan18_B5175_output.root",
	//	"HWPscan_jan19_B4353_output.root",
	//	"HWPscan_jan20_B3619_2_output.root",
	//	"HWPscan_jan20_B3619_output.root",
	//	"HWPscan_jan21_B3043_output.root",
	//	"HWPscan_jan22_B5173_output.root",
	//	"HWPscan_jan23_B5173_output.root"
	//};


	int Nfiles = filenames.size();

	TFile** f = new TFile*[Nfiles];

	TGraphErrors** g_blumlein = new TGraphErrors* [Nfiles];
	TGraphErrors** g_SNR = new TGraphErrors* [Nfiles];
	TGraphErrors** g_vibration = new TGraphErrors* [Nfiles];

	for (int i=0; i<Nfiles; i++){

		f[i] = TFile::Open(filenames[i]);
		g_blumlein[i] = (TGraphErrors*)f[i]->Get("HWPscan");
		g_SNR[i] = (TGraphErrors*)f[i]->Get("HWPscan_SNR");
		g_vibration[i] = (TGraphErrors*)f[i]->Get("HWPvibration");

		TString gTitle = filenames[i];
		gTitle.Remove(0,gTitle.Index("_jan")+1);
		gTitle.Remove(11);
		g_blumlein[i]->SetTitle(gTitle);
		g_SNR[i]->SetTitle(gTitle);
		g_vibration[i]->SetTitle(gTitle);
	}

	new TCanvas("","",1000,1000);
	for (int i=0; i<Nfiles; i++){
		g_blumlein[i]->GetXaxis()->SetLimits(0,60);
		g_blumlein[i]->GetYaxis()->SetRangeUser(-50,70);
		g_blumlein[i]->SetMarkerColor(i+1);
		g_blumlein[i]->Draw(i==0?"APL":"PL");
	}
	gPad->BuildLegend();


	new TCanvas("","",1000,1000);
	for (int i=0; i<Nfiles; i++){
		g_SNR[i]->GetXaxis()->SetLimits(0,60);
		g_SNR[i]->GetYaxis()->SetRangeUser(0,60);
		g_SNR[i]->SetMarkerColor(i+1);
		g_SNR[i]->Draw(i==0?"APL":"PL");
	}
	gPad->BuildLegend();

	new TCanvas("","",1000,1000);
	for (int i=0; i<Nfiles; i++){
		g_vibration[i]->GetXaxis()->SetLimits(0,60);
		g_vibration[i]->GetYaxis()->SetRangeUser(0,20);
		g_vibration[i]->SetMarkerColor(i+1);
		g_vibration[i]->Draw(i==0?"APL":"PL");
	}
	gPad->BuildLegend();














}
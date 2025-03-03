{


	TFile* f1 = TFile::Open("xyt_S12_noRF_unrandom_ext44_noVertCOD_each_calo_3Dhist.root");
	TFile* f2 = TFile::Open("xyt_S12_xRF_unrandom_ext44_noVertCOD_each_calo_3Dhist.root");
	TFile* f3 = TFile::Open("xyt_S12_xyRF_unrandom_ext_noVertCOD_each_calo_3Dhist.root");



	TH3F* h3_calo7_1 = (TH3F*)f1->Get("h3_xyt_EXT_calo7");
	h3_calo7_1->SetName("h3_calo7_1");
	TH3F* h3_calo7_2 = (TH3F*)f2->Get("h3_xyt_EXT_calo7");
	h3_calo7_2->SetName("h3_calo7_2");
	TH3F* h3_calo7_3 = (TH3F*)f3->Get("h3_xyt_EXT_calo7");
	h3_calo7_3->SetName("h3_calo7_3");

	TFile* fout = new TFile("beam_dists.root","recreate");

	TH2D* h2_calo7_1 = (TH2D*)h3_calo7_1->Project3D("yx");
	TH2D* h2_calo7_2 = (TH2D*)h3_calo7_2->Project3D("yx");
	TH2D* h2_calo7_3 = (TH2D*)h3_calo7_3->Project3D("yx");

	h2_calo7_1->SetName("noRF");
	h2_calo7_2->SetName("xRF");
	h2_calo7_3->SetName("xyRF");
	
	fout->Write();
	fout->Close();



}
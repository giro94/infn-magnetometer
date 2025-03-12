double bst = -45.;
double ben =  45.;
int bnum   = 180;

void makeBeamDistrib(TString fInput = "Beam_distribution_run3b.csv"){
    
    TTree *t = new TTree("t", "tree");
    t->ReadFile(fInput, "x:y:intensity");
    
    Float_t x, y, w;
    
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("intensity", &w);
    
    TH2F *h2_beam = new TH2F("h2_beam", "Beam Distribution; x[mm]; y[mm]", bnum, bst, ben, bnum, bst, ben);
    TH1F *h_beam_x = new TH1F("h_beam_x", "Beam Distribution; x[mm]; Intensity", bnum, bst, ben);
    TH1F *h_beam_y = new TH1F("h_beam_y", "Beam Distribution; y[mm]; Intensity", bnum, bst, ben);
    
    for (int i=0; i<t->GetEntries(); i++) {
        t->GetEntry(i);
        h2_beam->Fill(x, y, w);
        h_beam_x->Fill(x, w);
        h_beam_y->Fill(y, w);
    }
    
    
    double integral = h_beam_x->Integral(1, h_beam_x->GetNbinsX(), "width");
    TGraph *g_beam_x = new TGraph();
    
    for (int j=0; j<h_beam_x->GetNbinsX(); j++) {
        g_beam_x->SetPoint(j, h_beam_x->GetBinCenter(j+1), h_beam_x->GetBinContent(j+1)/integral);
    }
    
    
    TCanvas *c0 = new TCanvas("c0", "c0", 800, 600);
    h2_beam->Draw("colz");
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
    c1->Divide(2, 1);
    c1->cd(1);
    h_beam_x->Draw();
    c1->cd(2);
    h_beam_y->Draw();
    

    new TCanvas();
    g_beam_x->Draw();
    cout<<g_beam_x->Integral()<<endl;
    
    TFile *fOut = new TFile("BeamDistribution.root", "recreate");
    h2_beam->Write();
    h_beam_x->Write();
    h_beam_y->Write();
    g_beam_x->Write("g_beam_x");

    c1->Write();

    fOut->Write();
    fOut->Close();
}

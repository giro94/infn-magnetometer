{

    TFile* f_beam = TFile::Open("BeamDistrib/beam_dists.root");
    TFile* f_model = TFile::Open("KickerSpaceModel_1.root");


    TH2D* h2_beam = (TH2D*)f_beam->Get("noRF");
    //h2_beam->Scale(1./h2_beam->Interpolate(0,0));

    TH2D* h2_model = (TH2D*)f_model->Get("h2Space");

    for (int i=1; i<=h2_model->GetNbinsY(); i++) {
        for (int j=1; j<=h2_model->GetNbinsX(); j++) {
            double x=h2_model->GetXaxis()->GetBinCenter(j);
            double y=h2_model->GetYaxis()->GetBinCenter(i);
            if(x*x+y*y > pow(45., 2)) h2_model->SetBinContent(j, i, 0.);
        }
    }

    TH2D* h2_model_rebinned = (TH2D*)h2_beam->Clone("h2_model_rebinned");
    h2_model_rebinned->Reset();
    for (int i=1; i<=h2_model_rebinned->GetNbinsY(); i++) {
        for (int j=1; j<=h2_model_rebinned->GetNbinsX(); j++) {
            double x = h2_model_rebinned->GetXaxis()->GetBinCenter(j);
            double y = h2_model_rebinned->GetYaxis()->GetBinCenter(i);
            double z = h2_model->Interpolate(x,y);
            h2_model_rebinned->SetBinContent(j, i, z);
        }
    }
    
    TH2D* h2_convolution = (TH2D*)h2_beam->Clone("h2_convolution");
    h2_convolution->Reset();
    for (int i=1; i<=h2_convolution->GetNbinsX(); i++) {
        for (int j=1; j<=h2_convolution->GetNbinsY(); j++) {
            double beam = h2_beam->GetBinContent(i,j);
            double model = h2_model_rebinned->GetBinContent(i,j);
            h2_convolution->SetBinContent(i,j,beam*model);
        }
    }
    h2_convolution->Scale(1./h2_convolution->Interpolate(0,0));


    new TCanvas();
    h2_beam->Draw("colz");

    new TCanvas();
    h2_model->Draw("colz");

    new TCanvas();
    h2_model_rebinned->Draw("colz");

    new TCanvas();
    h2_convolution->Draw("colz");


}
{

    TFile* f_beam = TFile::Open("BeamDistrib/beam_dists_eva.root");
    TFile* f_model = TFile::Open("Tools/KickerSpaceModel_UMass.root");


    TH2D* h2_beam = (TH2D*)f_beam->Get("noRF");

    TH2D* h2_model = (TH2D*)f_model->Get("h2Space");

    for (int i=1; i<=h2_beam->GetNbinsY(); i++) {
        for (int j=1; j<=h2_beam->GetNbinsX(); j++) {
            double x=h2_beam->GetXaxis()->GetBinCenter(j);
            double y=h2_beam->GetYaxis()->GetBinCenter(i);
            if(x*x+y*y > pow(45., 2)) h2_beam->SetBinContent(j, i, 0.);
        }
    }
    h2_beam->Scale(1./h2_beam->Integral());

    for (int i=1; i<=h2_model->GetNbinsY(); i++) {
        for (int j=1; j<=h2_model->GetNbinsX(); j++) {
            double x=h2_model->GetXaxis()->GetBinCenter(j);
            double y=h2_model->GetYaxis()->GetBinCenter(i);
            if(x*x+y*y > pow(45., 2)) h2_model->SetBinContent(j, i, 0.);
        }
    }
    h2_model->Scale(1./h2_model->Interpolate(0,0));
    h2_model->Scale(1./0.782996);//0.741755);

    TH2D* h2_convolution = (TH2D*)h2_beam->Clone("h2_convolution");
    h2_convolution->Reset();
    for (int i=1; i<=h2_convolution->GetNbinsX(); i++) {
        for (int j=1; j<=h2_convolution->GetNbinsY(); j++) {
            double beam = h2_beam->GetBinContent(i,j);
            double model = h2_model->GetBinContent(i,j);
            h2_convolution->SetBinContent(i,j,beam*model);
        }
    }

    cout<<"Space factor: "<<h2_convolution->Integral()<<"\n";

    new TCanvas();
    h2_beam->Draw("colz");

    new TCanvas();
    h2_model->Draw("colz");

    new TCanvas();
    h2_convolution->Draw("colz");


}
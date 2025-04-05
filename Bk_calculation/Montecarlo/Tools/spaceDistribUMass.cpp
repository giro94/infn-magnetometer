{

TH2F *h2SpaceSiu = new TH2F("h2Space", "", 89, -110, 110, 65, -80, 80);

std::ifstream file("umass_2025_03_17.txt");
if (!file.is_open()) {
    std::cerr << "Errore: impossibile aprire il file " << std::endl;
    return;
}


std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double x, y, By;

        // Leggi le coordinate x e y dalla riga
        char delimiter;
        //if (ss >> x >> delimiter >> y >> delimiter >> Bx >> delimiter >> By && delimiter == '    ') {
        if (ss >> x >> y >> By) {
            h2SpaceSiu->Fill(x, y, By);
        } else {
            std::cerr << "Errore nel parsing della riga: " << line << std::endl;
        }
    }
h2SpaceSiu->Draw("colz");


TH2F *h2SpaceNorm = (TH2F*)h2SpaceSiu->Clone();
h2SpaceNorm->SetName("h2SpaceNorm"); 

double norm00 = -16.2 / 0.74; //-21.89
double norm10 = -35.4 / 0.678; //-52.21

double orig00 = h2SpaceSiu->Interpolate(0, 0);

h2SpaceNorm->Scale(1./orig00);

// Verifica dei valori normalizzati
cout << "Valori interpolati dopo la normalizzazione:\n"
     << "(0, 0): " << h2SpaceNorm->Interpolate(0, 0) << "\n"
     << "(17.5, 0): " << h2SpaceNorm->Interpolate(17.5, 0) << "\n"
     << "(-17.5, 0): " << h2SpaceNorm->Interpolate(-17.5, 0) << endl;

// Disegna il TH2 normalizzato
new TCanvas();
h2SpaceNorm->Draw("colz");


TFile *fOut = new TFile("KickerSpaceModel_UMass.root", "recreate");
h2SpaceNorm->Write("h2Space");
fOut->Write();
fOut->Close();


}
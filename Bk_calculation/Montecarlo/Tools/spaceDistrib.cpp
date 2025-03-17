int nSlices=21;


void spaceDistrib(){

    gStyle->SetOptStat(0000);

    TGraph *spx = new TGraph();
int i=0;

spx->SetPoint(i, -45.1, 0); i++;
spx->SetPoint(i, -42.73969179419197, 65.14270675169547); i++;
spx->SetPoint(i, -41.63201281902776, 66.39766170665047); i++;
spx->SetPoint(i, -40.5058432114496, 67.50480817896548); i++;
spx->SetPoint(i, -39.218792231360284, 68.5423625873064); i++;
spx->SetPoint(i, -37.529537819993046, 69.55461079056585); i++;
spx->SetPoint(i, -35.77995289393413, 70.5636957181901); i++;
spx->SetPoint(i, -33.970037453183517, 71.43675979350137); i++;
spx->SetPoint(i, -32.039460983049537, 71.6929851199514); i++;
spx->SetPoint(i, -30.108884512915557, 71.76890373519586); i++;
spx->SetPoint(i, -28.178308042781572, 71.44624962040692); i++;
spx->SetPoint(i, -26.247731572647592, 71.13308533252354); i++;
spx->SetPoint(i, -24.43781613189698, 70.46879744913453); i++;
spx->SetPoint(i, -22.56757017645469, 69.73175422613625); i++;
spx->SetPoint(i, -20.817985250395767, 68.74164895232312); i++;
spx->SetPoint(i, -19.00806980964516, 67.82113574248405); i++;
spx->SetPoint(i, -17.19815436889455, 66.93858184026723); i++;
spx->SetPoint(i, -15.388238928143947, 66.07500759186152); i++;
spx->SetPoint(i, -13.578323487393336, 65.23041299726692); i++;
spx->SetPoint(i, -11.768408046642724, 64.42377771029456); i++;
spx->SetPoint(i, -09.898162091200433, 63.649829604885774); i++;
spx->SetPoint(i, -08.088246650449822, 62.89907885413502); i++;
spx->SetPoint(i, -06.218000695007531, 62.41193440631642); i++;
spx->SetPoint(i, -04.2874242248735506, 62.165198906771934); i++;
spx->SetPoint(i, -02.3568477547395617, 62.0608108108108); i++;
spx->SetPoint(i, -00.4262712846055816, 62.01336167628301); i++;
spx->SetPoint(i, 01.5646357002200872, 62.01336167628301); i++;
spx->SetPoint(i, 03.555542685045756, 62.00387184937746); i++;
spx->SetPoint(i, 05.365458125796367, 62.35499544488307); i++;
spx->SetPoint(i, 07.296034595930347, 62.772547828727596); i++;
spx->SetPoint(i, 09.226611066064327, 63.3039781354388); i++;
spx->SetPoint(i, 11.096857021506619, 63.81010223706852); i++;
spx->SetPoint(i, 12.90677246225723, 64.72112562000201); i++;
spx->SetPoint(i, 14.656357388316152, 65.60051624658365); i++;
spx->SetPoint(i, 16.526603343758444, 66.43351216384923); i++;
spx->SetPoint(i, 18.336518784509055, 67.31922934170125); i++;
spx->SetPoint(i, 20.146434225259657, 68.1290279043088); i++;
spx->SetPoint(i, 21.89601915131858, 69.16869116307316); i++;
spx->SetPoint(i, 23.70593459206919, 70.07971454600667); i++;
spx->SetPoint(i, 25.515850032819802, 70.89583965988459); i++;
spx->SetPoint(i, 27.325765473570414, 71.56961737017915); i++;
spx->SetPoint(i, 29.256341943704385, 72.00614940783478); i++;
spx->SetPoint(i, 31.186918413838374, 72.07257819617368); i++;
spx->SetPoint(i, 33.177825398664034, 71.82373384620574); i++;
spx->SetPoint(i, 34.927410324722956, 71.00022775584571); i++;
spx->SetPoint(i, 36.737325765473567, 70.05124506528999); i++;
spx->SetPoint(i, 38.4265801768408, 68.86501670209535); i++;
spx->SetPoint(i, 40.05550407351635, 67.56385599082226); i++;
spx->SetPoint(i, 41.62409745550022, 66.15514390795288); i++;
spx->SetPoint(i, 43.13258096010324, 64.51311150636988); i++;
spx->SetPoint(i, 45.1, 0); i++;

double normPt = spx->Eval(0.);
for (int j = 0; j < spx->GetN(); ++j)
{
    spx->GetY()[j] *= -16.2/normPt;
}

spx->SetMarkerColor(kBlack);
spx->SetMarkerStyle(8);
spx->SetMarkerSize(0.5);
spx->Draw();
spx->GetYaxis()->SetRangeUser(-50, 0.1);

TF1 *fKick = new TF1("fKick", "[0] + [1]*x^2 + [2]*x^4", -45, 45);
//fKick->SetParameters(-16.1398  , -0.00039315, -0.00532216, 7.23506e-07, 2.70196e-06);
fKick->SetParameters(-16.1398, -0.00532216, 2.70196e-06);

spx->Fit(fKick, "RM");

cout<<spx->Eval(17.5)<<endl;

TGraph *g = new TGraph();
g->SetPoint(0, -45, 0);
g->SetPoint(1, -17.5, -35.4);
g->SetPoint(2, 0, -16.2);
g->SetPoint(3, 17.5, -35.4);
g->SetPoint(4, 45, 0);
//for (int j = 0; j < g->GetN(); ++j)
//{
//    g->GetY()[j] *= 1./16.2;
//}

g->SetMarkerColor(kRed);
g->SetMarkerStyle(8);
g->SetMarkerSize(0.5);

g->Draw("PSAME");

//g->Fit(fKick, "RM");

TF1 *fTr = new TF1("fTr", "[0] + [1]*x^2 + [2]*x^4", -45, 45);
fTr->SetParameters(-16.2, -0.13618820, 6.725342e-5);
fTr->SetLineColor(kBlue);
g->Fit(fTr, "RM0");
//fTr->Draw("same");
cout<<"ciao"<<endl;
cout<<fTr->Eval(17.5)<<endl;


TGraph *spx1 = (TGraph*)spx->Clone();

double normPt1 = spx1->Eval(0.);

for (int j = 0; j < spx1->GetN(); ++j)
{
    spx1->GetY()[j] = (spx1->GetY()[j]-normPt1) * (-35.4 - normPt1) /(spx->Eval(17.5)- normPt1) + normPt1;
}
spx1->SetMarkerColor(kMagenta);
spx1->Draw("SAME");

spx1->Fit(fTr, "RM");
fTr->Draw("same");

new TCanvas();

TH2F *h2SpaceSiu = new TH2F("h2Space", "", 89, -110, 110, 65, -80, 80);

std::ifstream file("KickerBfield2d_E989.csv");
if (!file.is_open()) {
    std::cerr << "Errore: impossibile aprire il file " << std::endl;
    return;
}


std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double x, y, Bx, By;

        // Leggi le coordinate x e y dalla riga
        char delimiter;
        //if (ss >> x >> delimiter >> y >> delimiter >> Bx >> delimiter >> By && delimiter == '    ') {
        if (ss >> x >> y >> Bx >> By) {
            h2SpaceSiu->Fill(x, y, By);
        } else {
            std::cerr << "Errore nel parsing della riga: " << line << std::endl;
        }
    }
h2SpaceSiu->Draw("colz");


TH2F *h2SpaceNorm = (TH2F*)h2SpaceSiu->Clone();
h2SpaceNorm->SetName("h2SpaceNorm"); 

double norm00 = -16.2/0.978; 
double norm10 = -35.4/0.993; 

double orig00 = h2SpaceSiu->Interpolate(0, 0);     
double orig10 = h2SpaceSiu->Interpolate(17.5, 0);  
double orig_10 = h2SpaceSiu->Interpolate(-17.5, 0);

if (std::abs(orig10 - orig_10) > 1e-6) {
    cerr << "Errore: i valori originali a (17.5, 0) e (-17.5, 0) non corrispondono." << endl;
    return;
}

for (int i = 1; i <= h2SpaceNorm->GetNbinsX(); ++i) { // Nota: gli indici partono da 1
    for (int j = 1; j <= h2SpaceNorm->GetNbinsY(); ++j) {
        double origVal = h2SpaceNorm->GetBinContent(i, j);

        double corrVal = (origVal - orig00) * (norm10 - norm00) / (orig10 - orig00) + norm00;

        h2SpaceNorm->SetBinContent(i, j, corrVal);
    }
}

// Verifica dei valori normalizzati
cout << "Valori interpolati dopo la normalizzazione:\n"
     << "(0, 0): " << h2SpaceNorm->Interpolate(0, 0) << "\n"
     << "(17.5, 0): " << h2SpaceNorm->Interpolate(17.5, 0) << "\n"
     << "(-17.5, 0): " << h2SpaceNorm->Interpolate(-17.5, 0) << endl;

// Disegna il TH2 normalizzato
new TCanvas();
h2SpaceNorm->Draw("colz");


TFile *fOut = new TFile("KickerSpaceModel_1.root", "recreate");
h2SpaceNorm->Write("h2Space");
fOut->Write();
fOut->Close();
     
}

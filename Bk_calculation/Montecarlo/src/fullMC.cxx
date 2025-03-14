#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <string>
#include <unistd.h>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TRandomGen.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

using namespace std;

#define NSLICESX 21
#define NSLICESY 21

int eInt(double e);
TH1D *gKick = nullptr;
TH1F *hIntegral = nullptr;
double B0 = 1.45e+7; //Magnet field in mG
double f_azimuth = 0.085; //Kickers cover only 8.5% of azimut
double f_kickers = (53.1+53.0+55.0)/(3*55.0); //average of kickers strength for K3
//double f_kickers = (53.1+53.0+55.0)/(3*53.1); //average of kickers strength for K1
double w0 = 1.43937;
double N0 = 5.e+10;
int nBins = 4692;
int storageRadius = 45.;

double fitStartTime=30.136;
double fitEndTime=650.1;
double spaceNormalizationFactor = 1.;

/*********************************************
 Usage:
 ./bin/fullMC <radius [0, 1]> <space model> <kicker transient hist name> <beam hist name>
 Radius 0 is at the magic radius, 1 is at 17.5 mm
 Space model can be:
    flat: flat space distribution, no kick weighting
    x2: x^2 distribution, y flat
    x4: x^4 distribution, y flat
    x4y2: x^4, y^2 distribution 
    x2m2: x^2 distrib with -2mm shift (for uncertainty)
    x2p2: x^2 distrib with +2mm shift (for uncertainty)
    x4m2: x^4 distrib with -2mm shift (for uncertainty)
    x4p2: x^4 distrib with +2mm shift (for uncertainty)
    hist: uses 2d histogram distribution in the file KickerSpaceModel_1.root named h2Space

 Kicker transient hist name follows the naming Paolo gave (see INFN_Umass_hd.root file):
    hist_name           Vibr    Smooth
    h1_kick1_R0_ra:     yes,    yes
    h1_kick1_R0_p:      no,     no
    h1_kick1_R0_p_ra:   no,     yes
    h1_kick1_R0:        yes,    no
    h1_kick1_R1_ra:     no,     yes
    h1_kick1_R1:        no,     no

 Beam hist name is the name of the beam distribution hist in the file BeamDistrib/beam_dists.root:
    noRF 
    xRF
    xyRF

 The output file is written following the parameters given:
 full_r<radius>_m<spacemodel>_K-<kickertransient>_B-<beamdist>.root 
 *********************************************/


Double_t fitFunc13Par(Double_t *x, Double_t *par){
    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    
    Double_t A = par[2];
    Double_t w = par[3];
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    Double_t ACBO = par[5];
    Double_t wCBO = par[6];
    while (par[7] < 0)              par[7] += TMath::TwoPi();
    while (par[7] > TMath::TwoPi()) par[7] -= TMath::TwoPi();
    Double_t pCBO = par[7];
    Double_t tCBO = par[8];
    Double_t AVW = par[9];
    Double_t wVW = par[10];
    //limit phase to be between 0 and 2pi
    while (par[11] < 0)              par[11] += TMath::TwoPi();
    while (par[11] > TMath::TwoPi()) par[11] -= TMath::TwoPi();
    Double_t pVW = par[11];
    Double_t kLM = par[12];

    
    Double_t CBOTerm = 1.0 + (exp(-t / tCBO) * ACBO * cos(wCBO * t + pCBO));
    Double_t tVW = tCBO * wCBO / wVW;
    Double_t VWTerm = 1.0 + exp(-t / tVW) * AVW * cos(wVW * t + pVW);
    
    Double_t LMTerm = 1.0 - kLM*(hIntegral->Interpolate(t));

    return N * exp(-t / B) * (1 + A * cos( w * t + p)) * CBOTerm * VWTerm * LMTerm;
}

Double_t fitFunc13Par_EC(Double_t *x, Double_t *par){
    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    
    Double_t A = par[2];
    Double_t w = par[3];
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    Double_t ACBO = par[5];
    Double_t wCBO = par[6];
    while (par[7] < 0)              par[7] += TMath::TwoPi();
    while (par[7] > TMath::TwoPi()) par[7] -= TMath::TwoPi();
    Double_t pCBO = par[7];
    Double_t tCBO = par[8];
    Double_t AVW = par[9];
    Double_t wVW = par[10];
    //limit phase to be between 0 and 2pi
    while (par[11] < 0)              par[11] += TMath::TwoPi();
    while (par[11] > TMath::TwoPi()) par[11] -= TMath::TwoPi();
    Double_t pVW = par[11];
    Double_t kLM = par[12];

    
    Double_t CBOTerm = 1.0 + (exp(-t / tCBO) * ACBO * cos(wCBO * t+ pCBO));
    Double_t tVW = tCBO * wCBO / wVW;
    Double_t VWTerm = 1.0 + exp(-t / tVW) * AVW * cos(wVW * t + pVW);
    
    Double_t LMTerm = 1.0 - kLM*(hIntegral->Interpolate(t));
    
    Double_t BkFactor = par[13];
    Double_t KickTerm = BkFactor * (gKick->Integral(0, gKick->FindBin(t), "width") * (f_azimuth * f_kickers) / B0); //convert mG to ppb
    
    return N * exp(-t / B) * (1 + A * cos( w * (t + KickTerm) + p)) * CBOTerm * VWTerm * LMTerm;
}

Double_t fitFunc5Par(Double_t *x, Double_t *par){
    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    
    Double_t A = par[2];
    Double_t w = par[3];
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];

    return N * exp(-t / B) * (1 + A * cos( w * t + p));
}

Double_t fitFunc5Par_EC(Double_t *x, Double_t *par){
    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    
    Double_t A = par[2];
    Double_t w = par[3];
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    
    Double_t BkFactor = par[5];
    Double_t KickTerm = BkFactor * (gKick->Integral(0, gKick->FindBin(t), "width") * (f_azimuth * f_kickers) / B0); //convert mG to ppb
    
    return N * exp(-t / B) * (1 + A * cos( w * (t + KickTerm) + p));
}

int main(int argc, char *argv[])
{
    int rs = stoi(argv[1]);
    string hn = argv[3];
    pair<double, double> normalizationPoint = {0., 0.};
    double spaceNormalizationFactor = 1.0;

    string sp=argv[2];
    if (sp=="") sp="flat";
    cout<<"Using "<<sp<<" spatial model"<<endl;
    
    TFile *intFile = TFile::Open("INFN_UMass_average_hd.root");
    gKick = (TH1D*)intFile->Get(hn.c_str());

    if(rs==0){
    //    TFile *intFile = TFile::Open("INFN_golden_R0_Bon_h.root");
    //    gKick = (TH1D*)intFile->Get("h_R0_Bon");
        normalizationPoint = {0., 0.};
        if(sp=="hist" or sp == "x4y2") spaceNormalizationFactor = 0.987;
    }
    else if (rs == 1){
    //    TFile *intFile = TFile::Open("INFN_golden_R1_Bon_h.root");
    //    gKick = (TH1D*)intFile->Get("h_R1_Bon");
        normalizationPoint = {17.5, 0.}; //usando 23. qui viene lo stesso risultato che per R0
        if(sp=="hist" or sp == "x4y2") spaceNormalizationFactor = 0.993;
    }

    
    bool use5par = false;
    int nPars = 13;
    int bkfactorPar = 13;
    
    if(use5par){
        nPars = 5;
        bkfactorPar = 5;
    }

    
    TFile *file = new TFile("fit_2D_eBinned_allFloat.root");
    TF1 *f1 = (TF1*)((TH1F*)file->Get("hWiggle_12par"))->GetFunction("f1");
    TFile *fLM = TFile::Open("Jt2D_Lt_nh1-1_dt125-125_x80-80_lkeff90.root");
    hIntegral = (TH1F*)fLM->Get("Jt");
    
    //string beamDistribFileName="BeamDistrib/BeamDistribution_3b.root";
    string beamDistribFileName="BeamDistrib/beam_dists.root";
    TFile *fBD = TFile::Open(beamDistribFileName.c_str());

    string beamHistName = "h2_beam";

    string bd = argc>4 ? argv[4]: "";
    if(bd != ""){
        beamHistName = bd;
    }
    cout<<"Using Beam Distribution: "<<beamHistName.c_str()<<endl;
    TH2D *h2Beam = (TH2D*) fBD->Get(beamHistName.c_str());
    cout<<"Integral of beam distribution: "<<h2Beam->Integral()<<endl;
    h2Beam->Scale(1./h2Beam->Integral());

    TH2D *h2Space = nullptr;

    if(sp=="hist"){
        TFile *fSP = TFile::Open("KickerSpaceModel_1.root");
        h2Space = (TH2D*) fSP->Get("h2Space");

        for (int i=1; i<=h2Space->GetNbinsY(); i++) {
            for (int j=1; j<=h2Space->GetNbinsX(); j++) {

                double x=h2Space->GetXaxis()->GetBinCenter(j);
                double y=h2Space->GetYaxis()->GetBinCenter(i);
                if(x*x+y*y > pow(storageRadius, 2)) h2Space->SetBinContent(j, i, 0.);
            }
        }
    }
    else{
        //La logica è di definire NSLICES funzioni in cui pesare N0 con la funzione di distribuzione del fascio, quindi N0/fSpace->Eval(r), mentre il fattore Bk con la parabola (il picco a 1 corrisponde a R0, quindi sarà normalizzato in questo modo)
        TF1 *fSpace_x = new TF1("fSpace_x", "[0] + [1]*x^2 + [2]*x^4", -45, 45); //x is radius in mm;
        TF1 *fSpace_y = new TF1("fSpace_y", "[0] + [1]*x^2", -45, 45); //x is radius in mm;
        fSpace_y->SetParameters(1., 0.);

        if(sp=="x4") fSpace_x->SetParameters(-16.2, -0.0752902, 4.1131e-05); //x^4  
        if(sp=="x4y2") {
            fSpace_x->SetParameters(-16.2, -0.0752902, 4.1131e-05); //x^4 
            fSpace_y->SetParameters(-1., -0.000273479);
        }
        if(sp=="flat"){
            fSpace_x->SetParameters(-1, 0, 0); 
            fSpace_y->SetParameters(-1., 0);
        }
        if(sp=="x2"){
            fSpace_x->SetParameters(-16.2, -0.063, 0); 
            fSpace_y->SetParameters(-1., 0);
        }
        if(sp=="x2all"){
            fSpace_x->SetParameters(-16.1, -0.0637, 0); 
            fSpace_y->SetParameters(-1., 0);
        }
        if(sp=="x2m2"){
            fSpace_x->SetParameters(-15.87, -0.0817, 0); 
            fSpace_y->SetParameters(-1., 0);
        }
        if(sp=="x2p2"){
            fSpace_x->SetParameters(-15.996, -0.0510, 0); 
            fSpace_y->SetParameters(-1., 0);
        }

        if(sp=="x4m2"){
            fSpace_x->SetParameters(-15.8249, -0.0939701, 5.02642e-5); 
            fSpace_y->SetParameters(-1., 0);
        }
        if(sp=="x4p2"){
            fSpace_x->SetParameters(-15.94, -0.0651522, 3.60611e-5); 
            fSpace_y->SetParameters(-1., 0);
        }

        
        h2Space = new TH2D("h2Space", "Bk XY Space Distribution", 920, -46, 46, 920, -46, 46);
        for (int i=1; i<=h2Space->GetNbinsY(); i++) {
            for (int j=1; j<=h2Space->GetNbinsX(); j++) {

                double x=h2Space->GetXaxis()->GetBinCenter(j);
                double y=h2Space->GetYaxis()->GetBinCenter(i);
                
                if(x*x+y*y <= pow(storageRadius, 2)) h2Space->SetBinContent(j, i, fSpace_x->Eval(x) * fSpace_y->Eval(y));
                else h2Space->SetBinContent(j, i, 0.);
            }
        }
    } 


    h2Space->Scale(1./h2Space->GetBinContent(h2Space->GetXaxis()->FindBin(normalizationPoint.first), h2Space->GetYaxis()->FindBin(normalizationPoint.second)));

    cout<<"Integral of space distribution: "<<h2Space->Integral()<<endl;
    cout<<"Value at 0: "<<h2Space->Interpolate(0., 0.)<<" and at Normalization point ("<<normalizationPoint.first<<", "<<normalizationPoint.second<<"): "<<h2Space->Interpolate(normalizationPoint.first, normalizationPoint.second)<<endl;


    TH1D *hSpace_x = (TH1D*) h2Space->ProjectionX("hSpace_x", h2Space->GetYaxis()->FindBin(0.), h2Space->GetYaxis()->FindBin(0.)); 
    TH1D *hSpace_y = (TH1D*) h2Space->ProjectionY("hSpace_y", h2Space->GetXaxis()->FindBin(0.), h2Space->GetXaxis()->FindBin(0.)); 


    map<int, map<int, TF1*>> fWiggle, fWiggle_EC;
    vector<double> radiusVec, yVec;
    
    double r_min = -storageRadius;
    double r_max = storageRadius;
    double delta_r = (r_max - r_min)/NSLICESX;
    
    double y_min = -storageRadius;
    double y_max = storageRadius;
    double delta_y = (y_max - y_min)/NSLICESY;

    double y = y_min+delta_y/2.;
    
    double N_total = 0;
    int bin_beam_x0, bin_beam_x1, bin_beam_y0, bin_beam_y1;
    int bin_space_x0, bin_space_x1, bin_space_y0, bin_space_y1;

    for (int i=0; i<NSLICESY; i++) {

        yVec.push_back(y);
        double radius = r_min+delta_r/2.;
        
        bin_beam_y0 = h2Beam->GetYaxis()->FindBin(y-delta_y/2.);
        bin_beam_y1 = h2Beam->GetYaxis()->FindBin(y+delta_y/2.) -1;
        bin_space_y0 = h2Space->GetYaxis()->FindBin(y-delta_y/2.);
        bin_space_y1 = h2Space->GetYaxis()->FindBin(y+delta_y/2.) -1;

        for(int j=0; j<NSLICESX; j++){

            if(j==0) radiusVec.push_back(radius);
            
            bin_beam_x0 = h2Beam->GetXaxis()->FindBin(radius-delta_r/2.);
            bin_beam_x1 = h2Beam->GetXaxis()->FindBin(radius+delta_r/2.) - 1;
            bin_space_x0 = h2Space->GetXaxis()->FindBin(radius-delta_r/2.);
            bin_space_x1 = h2Space->GetXaxis()->FindBin(radius+delta_r/2.) - 1;
            
            double N_r=N0 * h2Beam->Integral(bin_beam_x0, bin_beam_x1, bin_beam_y0, bin_beam_y1);
            double BkFactor = h2Space->Integral(bin_space_x0, bin_space_x1, bin_space_y0, bin_space_y1, "width") / (delta_r*delta_y);
            
            if(use5par){
                fWiggle_EC[i][j] = new TF1(Form("fWiggle_EC_%i_%i", i, j), fitFunc5Par_EC , 0., 700., 6); //parameter 14 is the Bk weight factor
                fWiggle[i][j] = new TF1(Form("fWiggle_%i_%i", i, j), fitFunc5Par , 0., 700., 5);
            }
            else{
                fWiggle_EC[i][j] = new TF1(Form("fWiggle_EC_%i_%i", i, j), fitFunc13Par_EC , 0., 700., 14); //parameter 14 is the Bk weight factor
                fWiggle[i][j] = new TF1(Form("fWiggle_%i_%i", i, j), fitFunc13Par , 0., 700., 13);
            }
            

            fWiggle_EC[i][j]->SetNpx(5000);
            fWiggle_EC[i][j]->SetParameter(3, w0);
            fWiggle_EC[i][j]->SetParameter(0, N_r);
            fWiggle_EC[i][j]->SetParameter(bkfactorPar, BkFactor);
            
            fWiggle[i][j]->SetNpx(5000);
            fWiggle[i][j]->SetParameter(3, w0);
            fWiggle[i][j]->SetParameter(0, N_r);
            
            for (int pn=1; pn<nPars; pn++) {
                if(pn==3) continue;
                fWiggle_EC[i][j]->SetParameter(pn, f1->GetParameter(pn));
                fWiggle[i][j]->SetParameter(pn, f1->GetParameter(pn));
            }

            
            
            //cout<<i<<" "<<y-delta_y/2.<< " "<<y+delta_y/2.<<" :: "<<j<<" "<<radius-delta_r/2.<< " "<<radius+delta_r/2.<<" :: "<<N_r/N0<<" "<<BkFactor<<endl;
            N_total+= N_r;
            
            radius += delta_r;
        }
        y += delta_y;
    }
    
    cout<<"Total number of positrons: "<<N_total<<endl;

    cout<<"Map Size: "<< fWiggle_EC.size()<<endl; 
    
    TString fOutName = "";
    if(bd == ""){
        fOutName = Form("full_r%s_m%s_K-%s.root", argv[1], argv[2], argv[3]);
    }
    else{
        fOutName = Form("full_r%s_m%s_K-%s_B-%s.root", argv[1], argv[2], argv[3], argv[4]);
    }

    TFile *fOut = new TFile(fOutName, "RECREATE");

    TH1F *hError_EC = new TH1F("hError_EC", "hError_EC", 1000, -50, 50);
    TH1F *hError = new TH1F("hError", "hError", 1000, -50, 50);
    TH1F *hDiff = new TH1F("hDiff", "hDiff", 1000, -50, 50);
    
    
    TH1F *hWiggle_EC = new TH1F("hWiggle_EC", "Wiggle Plot with Magnetometer Data; Time [#mus]; Counts [-]", nBins, 0., 700.);
    TH1F *hWiggle = new TH1F("hWiggle", "Wiggle Plot; Time [#mus]; Counts [-]", nBins, 0., 700.);
    
    hWiggle_EC->Sumw2();
    hWiggle->Sumw2();
    cout<<"Starting Loop\n";
    for (int k=1; k<=nBins; k++) {
        double t = hWiggle_EC->GetBinCenter(k);
        double N_EC = 0;
        double N = 0;

        for (int i=0; i<NSLICESY; i++) {
            for(int j=0; j<NSLICESX; j++){
                //N_beam = h3Beam->GetBinContent(x, y, t);
                N_EC += fWiggle_EC[i][j]->Eval(t); //* N_beam;
                N += fWiggle[i][j]->Eval(t);
            }
        }

        double dx = N; //rand.PoissonD(N);
        hWiggle->SetBinContent(k, dx);
        hWiggle->SetBinError(k, sqrt(dx));
        double dy = N_EC; //rand.PoissonD(N_EC);
        hWiggle_EC->SetBinContent(k, dy);
        hWiggle_EC->SetBinError(k, sqrt(dy));
    }

    
    //set the fitting function
    TF1 *fWiggle_fit = nullptr;
    
    if(use5par){
        fWiggle_fit = new TF1("fWiggle_fit", fitFunc5Par , 30., 651., 5);
    }
    else{
        fWiggle_fit = new TF1("fWiggle_fit", fitFunc13Par , 30., 651., 13);
    }
    
    
    fWiggle_fit->SetNpx(5000);
    
    fWiggle_fit->SetParameter(0, N_total);
    fWiggle_fit->SetParameter(3, w0);
    for (int pn=1; pn<nPars; pn++) {
        if(pn==3) continue;
        fWiggle_fit->SetParameter(pn, f1->GetParameter(pn));
    }

    
    cout<<"Fitting Bk altered histogram\n";
    hWiggle_EC->Fit(fWiggle_fit, "M", "", fitStartTime, fitEndTime);
    double wfit_EC = fWiggle_fit->GetParameter(3);
    double error_EC =(wfit_EC - w0)/w0*1e+9;

    //reset again
    cout<<"Fitting normal histogram\n";
    fWiggle_fit->SetParameter(0, N_total);
    fWiggle_fit->SetParameter(3, w0);
    for (int pn=1; pn<nPars; pn++) {
        if(pn==3) continue;
        fWiggle_fit->SetParameter(pn, f1->GetParameter(pn));
    }
    
    hWiggle->Fit(fWiggle_fit, "M", "", fitStartTime, fitEndTime);
    double wfit = fWiggle_fit->GetParameter(3);
    double error = (wfit - w0)/w0*1e+9;
    
    hError->Fill(error);
    hError_EC->Fill(error_EC);
    hDiff->Fill(error_EC-error);
    
    cout<<error<<" "<<error_EC<<endl;

    h2Space->Write();
    h2Beam->Write();
    //fSpace_x->Write();
    //fSpace_y->Write();
    hSpace_x->Write();
    hSpace_y->Write();
    gKick->Write();

    fOut->Write();
    fOut->Close();
    return 0;
}

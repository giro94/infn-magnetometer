#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TKey.h>
#include <iostream>
#include <map>
#include <string>

void dynConvertTH1(string inputFileName="INFN.root") {
    // Apri il file di input
    TFile *inputFile = TFile::Open(inputFileName.c_str());
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Errore nell'apertura del file di input!" << std::endl;
        return;
    }

    // Ottieni la lista degli oggetti nel file
    TList *listOfKeys = inputFile->GetListOfKeys();
    if (!listOfKeys) {
        std::cerr << "Nessun oggetto trovato nel file di input!" << std::endl;
        inputFile->Close();
        return;
    }

    // Mappa per gli istogrammi di input
    std::map<int, TH1D*> inputHistograms;
    int index = 0;
    for (TObject *obj : *listOfKeys) {
        TKey *key = dynamic_cast<TKey*>(obj);
        if (!key) continue;

        std::cout << "Nome da TKey: " << key->GetName() << std::endl;

        TObject *genericObj = key->ReadObj();
        if (!genericObj) {
            std::cerr << "Errore: impossibile leggere l'oggetto \"" << key->GetName() << "\"." << std::endl;
            continue;
        }

        std::cout << "Classe dell'oggetto letto: " << genericObj->ClassName() << std::endl;

        if (genericObj->InheritsFrom(TH1::Class())) {
            TH1D *hist = dynamic_cast<TH1D*>(genericObj);
            if (hist) {
                // Forza il nome dell'istogramma a essere quello della chiave
                hist->SetName(key->GetName());
                inputHistograms[index++] = hist;
                std::cout << "Caricato istogramma: " << hist->GetName() << std::endl;
            } else {
                std::cerr << "Errore: L'oggetto \"" << key->GetName()
                << "\" non è un TH1D valido." << std::endl;
            }
        } else {
            std::cerr << "Oggetto \"" << key->GetName() << "\" non è un istogramma!" << std::endl;
        }
    }

    // Crea il file di output
    string outputFileName = inputFileName.substr(0, inputFileName.find_last_of('.')) + "_hd.root";

    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Errore nella creazione del file di output!" << std::endl;
        inputFile->Close();
        return;
    }

    // Mappa per i nuovi istogrammi
    std::map<int, TH1D*> outputHistograms;

    // Crea i nuovi istogrammi corrispondenti
    for (const auto &[idx, hist] : inputHistograms) {
        TH1D *newHist = new TH1D(hist->GetName(), Form("%s;Time [#mus]; Bk [mG]", hist->GetName()), 4962, 0, 700);
        outputHistograms[idx] = newHist;
    }

    // Esegui l'interpolazione e riempi i nuovi istogrammi
    for (const auto &[idx, newHist] : outputHistograms) {
        TH1D *oldHist = inputHistograms[idx];
        for (int bx = 1; bx <= newHist->GetNbinsX(); ++bx) {
            double x = newHist->GetBinCenter(bx);
            double val = oldHist->Interpolate(1e-3 * x);
            if (x < 30) val = 0;
            newHist->SetBinContent(bx, val);
        }
        
        // Disegna l'istogramma
        TCanvas *canvas = new TCanvas();
        newHist->Draw();
    }

    // Scrivi i nuovi istogrammi sul file di output
    outputFile->Write();
    outputFile->Close();
    inputFile->Close();
}

#include "../analysis_tools.C"

void plot_compare_QWP(){

	vector<TString> filenames = {
		"analysis/analysis_EC_jan23_B5173_H25_Q00.root",
		"analysis/analysis_EC_jan24_B5173_H25_Q62.5.root",
		"analysis/analysis_EC_jan25_B5173_H25_Q19.5.root",
		"analysis/analysis_EC_jan25_B5173_H25_Q22.5.root",
		"analysis/analysis_EC_jan26_B5173_H25_Q130.root",
		"analysis/analysis_EC_jan27_B5173_H25_Q85.root",
		"analysis/analysis_EC_jan28_B5173_H25_Q00.root"
	};

	int Nfiles = filenames.size();

	TFile** f = new TFile*[Nfiles];

	TH1D** h1_traces = new TH1D* [Nfiles];
	TH1D** h1_traces_ra = new TH1D* [Nfiles];

	for (int i=0; i<Nfiles; i++){

		f[i] = TFile::Open(filenames[i]);
		h1_traces[i] = ((TProfile*)f[i]->Get("trace"))->ProjectionX();
		h1_traces_ra[i] = (TH1D*)f[i]->Get("trace_ra");

		TString histTitle = filenames[i];
		histTitle.Remove(0,histTitle.Index("EC_")+3);
		h1_traces[i]->SetTitle(histTitle);
		h1_traces_ra[i]->SetTitle(histTitle+" (smooth)");
		h1_traces[i]->SetLineWidth(2);
		h1_traces_ra[i]->SetLineWidth(2);
	}


	gStyle->SetOptStat(0);
	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		h1_traces[i]->SetLineColor(i+1);
		h1_traces[i]->Draw("HIST SAME");
	}
	gPad->BuildLegend();

	new TCanvas();
	TH1D* h1_subtract = (TH1D*)h1_traces[0]->Clone("subtract");
	h1_subtract->Add(h1_traces[1],-1);
	h1_subtract->Scale(0.5);
	h1_subtract->Draw("HIST");


	new TCanvas();
	TH1D* h1_Q00_Q130 = (TH1D*)h1_traces_ra[4]->Clone("h1_Q00_Q130");
	h1_Q00_Q130->Add(h1_traces_ra[6],-1);
	h1_Q00_Q130->SetTitle("Q130 - Q00");
	h1_Q00_Q130->Draw("HIST");
	h1_traces_ra[4]->Draw("HIST SAME");


	new TCanvas();
	TH1D* h1_Q00_Q22p5 = (TH1D*)h1_traces[6]->Clone("h1_Q00_Q22p5");
	h1_Q00_Q22p5->Add(h1_traces[3],-1);
	h1_Q00_Q22p5 = runningAverage_5_10_15(h1_Q00_Q22p5);
	h1_Q00_Q22p5->SetTitle("Q00 - Q22.5");
	h1_Q00_Q22p5->Draw("HIST");
	h1_traces_ra[3]->Draw("HIST SAME");
	h1_traces_ra[6]->Draw("HIST SAME");

	gStyle->SetPalette(kRainBow);
	//Show the first 8 kicks
	TH1D** h1_kicks = new TH1D* [8];
	for (int i=0; i<8; i++){
		h1_kicks[i] = (TH1D*)f[0]->Get(Form("trace_kick%d",i+1));
	}
	TCanvas* can_kicks = new TCanvas("can_kicks","",1200,1200);
	can_kicks->Divide(1,3);
	for (int i=0; i<8; i++){
		h1_kicks[i]->SetLineWidth(i==0?2:1);
		h1_kicks[i]->SetLineColor(i==0?kRed:kBlue);
		h1_kicks[i]->GetXaxis()->SetTitleSize(0.05);
		h1_kicks[i]->GetXaxis()->SetLabelSize(0.05);
		
		can_kicks->cd(1);
		h1_kicks[i]->GetXaxis()->SetRangeUser(-0.6,0.4);
		h1_kicks[i]->GetYaxis()->SetRangeUser(-80,40);
		h1_kicks[i]->DrawCopy(TString("HIST")+(i==0?"":" SAME"));
		gPad->SetGridy();
		if(i==7)gPad->BuildLegend(0.91,0.1,0.995,0.9);

		can_kicks->cd(2);
		h1_kicks[i]->GetXaxis()->SetRangeUser(0,0.7);
		h1_kicks[i]->GetYaxis()->SetRangeUser(-40,40);
		h1_kicks[i]->DrawCopy(TString("HIST")+(i==0?"":" SAME"));
		gPad->SetGridy();

		can_kicks->cd(3);
		h1_kicks[i]->GetXaxis()->SetRangeUser(3,4);
		h1_kicks[i]->GetYaxis()->SetRangeUser(-50,50);
		h1_kicks[i]->DrawCopy(TString("HIST")+(i==0?"":" SAME"));
		gPad->SetGridy();
	}



	//Show files 1 and 7 similarity Q00 Q00
	new TCanvas("","",1000,600);
	h1_traces[0]->SetLineColor(kBlack);
	h1_traces[6]->SetLineColor(kRed);
	h1_traces[0]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[6]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[0]->GetYaxis()->SetRangeUser(-60,40);
	h1_traces[6]->GetYaxis()->SetRangeUser(-60,40);
	h1_traces[0]->DrawCopy("HIST");
	h1_traces[6]->DrawCopy("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.4,0.2,0.85,0.35);



	//Show files 2 and 5 similarity Q62.5 Q130
	new TCanvas("","",1000,600);
	h1_traces[1]->SetLineColor(kBlack);
	h1_traces[4]->SetLineColor(kRed);
	h1_traces[1]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[4]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[1]->GetYaxis()->SetRangeUser(-40,70);
	h1_traces[4]->GetYaxis()->SetRangeUser(-40,70);
	h1_traces[1]->DrawCopy("HIST");
	h1_traces[4]->DrawCopy("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.4,0.7,0.85,0.85);



	//Show Q00 Q22.5 Q130
	new TCanvas("","",1000,600);
	h1_traces[3]->SetLineColor(kBlack);
	h1_traces[4]->SetLineColor(kRed);
	h1_traces[6]->SetLineColor(kBlue);
	h1_traces[3]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[4]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[6]->GetXaxis()->SetRangeUser(4.5,6);
	h1_traces[3]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces[4]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces[6]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces[3]->DrawCopy("HIST");
	h1_traces[4]->DrawCopy("HIST SAME");
	h1_traces[6]->DrawCopy("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.4,0.7,0.85,0.85);


	//Show Q00 Q22.5 Q130
	new TCanvas("","",1000,600);
	h1_traces_ra[3]->SetLineColor(kBlack);
	h1_traces_ra[4]->SetLineColor(kRed);
	h1_traces_ra[6]->SetLineColor(kBlue);
	h1_traces_ra[3]->GetXaxis()->SetRangeUser(4.5,8);
	h1_traces_ra[4]->GetXaxis()->SetRangeUser(4.5,8);
	h1_traces_ra[6]->GetXaxis()->SetRangeUser(4.5,8);
	h1_traces_ra[3]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces_ra[4]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces_ra[6]->GetYaxis()->SetRangeUser(-70,70);
	h1_traces_ra[3]->DrawCopy("HIST");
	h1_traces_ra[4]->DrawCopy("HIST SAME");
	h1_traces_ra[6]->DrawCopy("HIST SAME");
	gPad->SetGridy();
	gPad->BuildLegend(0.4,0.7,0.85,0.85);


	//Show Q00 Q22.5 Q130 zoomed window
	new TCanvas("","",1000,600);
	h1_traces_ra[3]->SetLineColor(kBlack);
	h1_traces_ra[4]->SetLineColor(kRed);
	h1_traces_ra[6]->SetLineColor(kBlue);
	h1_traces_ra[3]->GetXaxis()->SetRangeUser(6.4,6.8);
	h1_traces_ra[4]->GetXaxis()->SetRangeUser(6.4,6.8);
	h1_traces_ra[6]->GetXaxis()->SetRangeUser(6.4,6.8);
	h1_traces_ra[3]->GetYaxis()->SetRangeUser(-20,20);
	h1_traces_ra[4]->GetYaxis()->SetRangeUser(-20,20);
	h1_traces_ra[6]->GetYaxis()->SetRangeUser(-20,20);
	h1_traces_ra[3]->DrawCopy("HIST");
	h1_traces_ra[4]->DrawCopy("HIST SAME");
	h1_traces_ra[6]->DrawCopy("HIST SAME");
	gPad->SetGridy();
	//gPad->BuildLegend(0.4,0.7,0.85,0.85);



	for (int i=0; i<Nfiles; i++){
		new TCanvas();
		h1_traces[i]->Draw("HIST");
		gPad->BuildLegend();
	}



	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		for (int bx=1; bx<=h1_traces[i]->GetNbinsX(); bx++){
			h1_traces[i]->SetBinContent(bx,h1_traces[i]->GetBinContent(bx)-i*20);
		}
		h1_traces[i]->SetLineColor(i+1);
		h1_traces[i]->Draw("HIST SAME");
	}
	gPad->BuildLegend();


	new TCanvas();
	for (int i=0; i<Nfiles; i++){
		for (int bx=1; bx<=h1_traces_ra[i]->GetNbinsX(); bx++){
			h1_traces_ra[i]->SetBinContent(bx,h1_traces_ra[i]->GetBinContent(bx)-i*20);
		}
		h1_traces_ra[i]->SetLineColor(i+1);
		h1_traces_ra[i]->Draw("HIST SAME");
	}
	gPad->BuildLegend();

}
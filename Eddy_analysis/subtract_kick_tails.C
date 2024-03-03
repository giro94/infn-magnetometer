#include "../analysis_tools.C"

void subtract_kick_tails(TString input_file, TString output_file=""){

	TFile* f = TFile::Open(input_file);


	bool saveOutput = false;
	if (output_file != ""){
		saveOutput = true;
	}


	TH1D* trace = ((TProfile*)f->Get("trace_SNR2_timealigned"))->ProjectionX();
	TH1D** kicks = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks[i] = ((TProfile*)f->Get(Form("trace_SNR2_kick%d",i+1)))->ProjectionX();
	}
	TH1D* kick8long = ((TProfile*)f->Get("trace_SNR2_kick8long"))->ProjectionX();


	double kick_trigger = -200;
	vector<double> kick_timings;
	vector<double> kick_deltas;
	for (int bn=2; bn<=trace->GetNbinsX(); bn++){
		double time1 = trace->GetBinCenter(bn-1);
		double time2 = trace->GetBinCenter(bn);
		double val1 = trace->GetBinContent(bn-1);
		double val2 = trace->GetBinContent(bn);
		if (val1 > kick_trigger && val2 < kick_trigger){
			double fraction = (kick_trigger-val1)/(val2-val1);
			double kick_time = time1 + (time2-time1)*fraction;
			double kick_dt = kick_timings.size()<1 ? 0 : kick_time - kick_timings.back();
			kick_timings.push_back(kick_time);
			kick_deltas.push_back(kick_dt);
			cout<<"Found kick "<<kick_timings.size()<<" at "<<kick_timings.back()<<" ms (dt: "<<kick_deltas.back()<<" ms)\n";
		}
	}
	if (kick_timings.size() > 8){
		cout<<"ERROR: found more than 8 kicks\n";
		return;
	}


	for (int bn=1; bn<=trace->GetNbinsX(); bn++){
		if (trace->GetBinContent(bn) < -50){
			trace->SetBinContent(bn,0);
		}
	}

	for (int i=0; i<8; i++){
		for (int bn=1; bn<=kicks[i]->GetNbinsX(); bn++){
			if (kicks[i]->GetBinContent(bn) < -50){
				kicks[i]->SetBinContent(bn,0);
			}
		}
	}

	for (int bn=1; bn<=kick8long->GetNbinsX(); bn++){
		if (kick8long->GetBinContent(bn) < -50){
			kick8long->SetBinContent(bn,0);
		}
	}


	TH1D* trace_subtract = (TH1D*)trace->Clone("trace_subtract");
	TH1D** kicks_subtract = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_subtract[i] = (TH1D*)kicks[i]->Clone(Form("%s_subtract",kicks[i]->GetName()));
	}

	TH1D* trace_ra = runningAverage(runningAverage(runningAverage(trace,5,true),10,true),15,true);
	trace_ra->SetTitle(Form("%s (RunningAvg %d)",trace->GetTitle(),51015));	
	TH1D** kicks_ra = new TH1D*[8];
	for (int i=0; i<8; i++){
		kicks_ra[i] = runningAverage(runningAverage(runningAverage(kicks[i],5,true),10,true),15,true);
		kicks_ra[i]->SetTitle(Form("%s (RunningAvg %d)",kicks[i]->GetTitle(),51015));	
	}
	TH1D* kick8long_ra = runningAverage(runningAverage(runningAverage(kick8long,5,true),10,true),15,true);
	kick8long_ra->SetTitle(Form("%s (RunningAvg %d)",kick8long->GetTitle(),51015));	



	new TCanvas();
	trace->SetLineWidth(2);
	trace_ra->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_ra->SetLineColor(kRed);
	trace->Draw("HIST");
	trace_ra->Draw("HIST SAME");


	//Subtract fulltrace
	for (int i=7; i>=1; i--){
		//Subtract kick i-1 tail from kick i
		for (int bn=1; bn<=trace->GetNbinsX(); bn++){
			double this_time = trace->GetBinCenter(bn);

			if (this_time > kick_timings[i]){ // kick 8
				double tail_time = kick_timings[7] + this_time - kick_timings[i-1];
				double tail_bin = trace->FindBin(tail_time);
				if (tail_bin > trace->GetNbinsX()) break;
				if (tail_time > kick_timings[7] + 2*kick_deltas[7]) break;

				double tail_val = trace->GetBinContent(tail_bin);
				double this_val = trace_subtract->GetBinContent(bn);
				trace_subtract->SetBinContent(bn, this_val - tail_val);
			}
		}
	}

	//Cutout kicks
	for (int i=0; i<8; i++){
		for (int bn=1; bn<=kicks_subtract[i]->GetNbinsX(); bn++){
			double this_time = kicks_subtract[i]->GetBinCenter(bn);
			double trace_time = kick_timings[i]+this_time;
			int trace_bin = trace_subtract->FindBin(trace_time);
			double this_val = trace_subtract->GetBinContent(trace_bin);
			kicks_subtract[i]->SetBinContent(bn,this_val);
		}
	}


	new TCanvas();
	kicks[6]->SetLineWidth(2);
	kicks[7]->SetLineWidth(2);
	kicks[6]->SetLineColor(kBlue);
	kicks[7]->SetLineColor(kRed);
	kicks[6]->Draw("HIST");
	kicks[7]->Draw("HIST SAME");

	new TCanvas();
	trace->SetLineWidth(2);
	trace_subtract->SetLineWidth(2);
	trace->SetLineColor(kBlue);
	trace_subtract->SetLineColor(kRed);
	trace->Draw("HIST");
	trace_subtract->Draw("HIST SAME");


	TCanvas* can_kick_subtract = new TCanvas();
	can_kick_subtract->Divide(2,1);
	can_kick_subtract->cd(1);
	for (int i=7; i>=0; i--){
		kicks[i]->GetXaxis()->SetRangeUser(0,0.1);
		kicks[i]->GetYaxis()->SetRangeUser(-40,10);
		kicks[i]->SetLineWidth(2);
		kicks[i]->SetLineColor(i==0?kBlack:kRed);
		kicks[i]->Draw("HIST SAME");
	}
	can_kick_subtract->cd(2);
	for (int i=7; i>=0; i--){
		kicks_subtract[i]->GetXaxis()->SetRangeUser(0,0.1);
		kicks_subtract[i]->GetYaxis()->SetRangeUser(-40,10);
		kicks_subtract[i]->SetLineWidth(2);
		kicks_subtract[i]->SetLineColor(i==0?kBlack:kRed);
		kicks_subtract[i]->Draw("HIST SAME");
	}


}
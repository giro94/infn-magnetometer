#include <dirent.h>
#include "../analysis_tools.C"



void plot_aligned_kicks(TString folder = "HWP-0p5_QWPno_kicker667_nofilter_field_26Aug2022", TString output_folder = "./", int Nfilesmax = -1){

	bool SaveFigs = true;

	int Npoints= 195315; 
	float tstart = 0.;
	float tend = 100.;
	float dt = (tend-tstart)/Npoints;

	float t_before = -0.5;
	float t_after = 8.0;
	int Nkickbins = (t_after-t_before)/dt;
	float t_after8 = 22.0;
	int Nkickbins8 = (t_after8-t_before)/dt;

	double first_kick_guess = 4.5;
	double kick_dt_guess = 10.0;
	double kick_trigger = -200.0;


	TFile* fout = new TFile(Form("%s/output.root",output_folder.Data()),"recreate");

	TProfile* g_trace = new TProfile("trace","Trace;Time [ms];Voltage [mV]",Npoints,0,100);
	TProfile** g_trace_kicks = new TProfile*[8];
	for (int i=0; i<8; i++){
		TString hname = Form("trace_kick%d",i+1);
		TString htitle = Form("Trace Kick %d;Time [ms];Voltage [mV]",i+1);
		if (i==7){
			g_trace_kicks[i] = new TProfile(hname,htitle,Nkickbins8,t_before,t_after8);
		} else {
			g_trace_kicks[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
		}
	}
	TProfile* g_trace_avg = (TProfile*)g_trace_kicks[0]->Clone("trace_avg");
	g_trace_avg->Reset();
	g_trace_avg->SetTitle("Kick average");

	DIR* dir = opendir(folder.Data());
	struct dirent* dirfile;
	//for (int fi=2; fi<=254; fi++){
	//for (int fi=2; fi<=2; fi++){
	int fi=0;
	cout<<"Reading folder "<<folder<<"\n";

	vector<TString> files;
	while((dirfile = readdir(dir)) != NULL){
		if (dirfile->d_type != DT_REG) continue;
		TString fname = dirfile->d_name;
		if (!fname.EndsWith("csv")) continue;
		files.push_back(fname);
	}

	sort(files.begin(),files.end());

	int Nfiles = files.size();
	for (int fi=0; fi<Nfiles; fi++){

		if (Nfilesmax > 0 && fi >= Nfilesmax) break;

		TString fname = files[fi];
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";
		vector<pair<double,double>> trace;

		ifstream file;
		TString filename = Form("%s/%s",folder.Data(),fname.Data());
		file.open(filename.Data());
		if (!file.is_open()){
			cout<<"cannot open "<<filename<<"\n";
			return;
		}
		char buff[256];
		file.getline(buff,256); //Header, names

		stringstream ss(buff);
		int Nvar = 0;
		double Cpos = -1;
		while(ss.good()){
			string substr;
			getline(ss,substr,',');
			//cout<<Nvar<<" : "<<substr<<"\n";
			if (substr.find("average(C)") != string::npos){
				Cpos = Nvar;
			}
			Nvar++;
		}

		file.getline(buff,256); //Header, units
		file.getline(buff,256); //Header, blank row

		double time, A, B, C, avgC;
		char comma;
		while (!file.eof()){
			double var=0;
			for (int i=0; i<Nvar; i++){
				if (i>0) file>>comma;
				file>>var;

				if (i==0){time = var;}
				if (i==Cpos){avgC = var;}
			}

			if (file.eof()) break;

			trace.push_back(make_pair(time,avgC));
		}
		file.close();

		for (auto p : trace){
			time = p.first;
			avgC = p.second;
			g_trace->Fill(time,avgC);
		}



		for (int i=0; i<8; i++){
			double this_guess = first_kick_guess+i*kick_dt_guess;

			for (int j=0; j<trace.size(); j++){
				double time = trace[j].first;
				double avgC = trace[j].second;
				if (time > this_guess){
					double prev_time = trace[j-1].first;
					double prev_avgC = trace[j-1].second;
					if (prev_avgC > kick_trigger){
						if (avgC < kick_trigger){
							double fraction = (kick_trigger-prev_avgC)/(avgC-prev_avgC);
							double kick_time = prev_time + dt*fraction;

							//cout<<"Found kick "<<i+1<<" for trace "<<fi<<". ";
							//cout<<"Kick time: "<<kick_time<<", previous point: "<<prev_avgC<<"\n";

							for (int k=0; k<trace.size(); k++){
								double t_time = trace[k].first;
								double t_avgC = trace[k].second;
								if (t_time >= kick_time + t_before && t_time < kick_time + (i==7?t_after8:t_after)){
									g_trace_kicks[i]->Fill(t_time-kick_time,t_avgC);
								}
							}

							break;
						}
					}
				}
			}
		}
		fi++;
	}


	for (int i=0; i<8; i++){
		g_trace_avg->Add(g_trace_kicks[i]);
	}

	TCanvas* can_kicks = new TCanvas("can_kicks","canvas",1920,1080);
	can_kicks->Divide(3,3);

	for (int i=0; i<8; i++){
		can_kicks->cd(i+1);
		g_trace_kicks[i]->GetXaxis()->SetRangeUser(0,1);
		g_trace_kicks[i]->GetYaxis()->SetRangeUser(-20,10);
		g_trace_kicks[i]->Draw("HIST L");
	}
	can_kicks->cd(9);
	g_trace_avg->GetXaxis()->SetRangeUser(0,1);
	g_trace_avg->GetYaxis()->SetRangeUser(-20,10);
	g_trace_avg->Draw("HIST L");
	if(SaveFigs)can_kicks->SaveAs(Form("%s/all_kicks.png",output_folder.Data()));

	new TCanvas("","",1920,1080);
	g_trace_avg->GetXaxis()->SetRangeUser(0,1);
	g_trace_avg->Draw("HIST L");
	if(SaveFigs)gPad->SaveAs(Form("%s/averagekick.png",output_folder.Data()));

	new TCanvas("","",1920,1080);
	//g_trace->Rebin(8);
	g_trace->GetYaxis()->SetRangeUser(-80,80);
	g_trace->Draw("HIST L");
	if(SaveFigs)gPad->SaveAs(Form("%s/8kicks.png",output_folder.Data()));


	double fft_xmin = 0.1; //ms from kick
	double fft_xmax = 1.0; 
	new TCanvas("","",1920,1080);
    TH1 *fft_histogram = 0;
    TVirtualFFT::SetTransform(0);
    TH1F* fftResidualInit = SetupFFT(g_trace_avg, fft_xmin, fft_xmax);
    fft_histogram = fftResidualInit->FFT(fft_histogram,"MAG");
    TH1F* fftResidual = RescaleAxis(fft_histogram, 1./(fft_xmax - fft_xmin));
    fftResidual->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
    fftResidual->SetStats(0);
    fftResidual->SetName(Form("residualFFT [%.2f,%.2f] ms",fft_xmin, fft_xmax));
    fftResidual->Scale(1.0 / fftResidual->Integral());
    fftResidual->GetXaxis()->SetRangeUser(0, fftResidual->GetXaxis()->GetXmax()/2.);
    fftResidual->Draw("HIST");
    fftResidual->GetXaxis()->SetRangeUser(2,3000);
    gPad->SetLogx();
    fftResidual->GetXaxis()->SetRangeUser(2,3000);
	if(SaveFigs)gPad->SaveAs(Form("%s/FFT.png",output_folder.Data()));


	fout->Write();
	fout->Close();
}

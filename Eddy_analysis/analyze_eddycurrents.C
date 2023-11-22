#include <dirent.h>

vector<TString> getListOfFiles(TString folder){
	vector<TString> files;

	cout<<"Reading folder "<<folder<<"\n";
	DIR* dir = opendir(folder.Data());
	if (!dir) {
		cout<<"Cannot open folder "<<folder<<"!\n";
		throw std::runtime_error{"Cannot open folder"};
	}

	struct dirent* dirfile;
	while((dirfile = readdir(dir)) != NULL){
		if (dirfile->d_type != DT_REG) continue;
		TString fname = dirfile->d_name;
		if (!fname.EndsWith("csv")) continue;
		files.push_back(fname);
	}
	sort(files.begin(),files.end());

	return files;
}

TDatime getFileTime(TString filename){
	filename.ReplaceAll("SD_","");
	filename.ReplaceAll(".csv","");
	std::tm tm{};
	std::istringstream iss(filename.Data());
	iss >> std::get_time(&tm, "%m_%d_%Y %H_%M_%S");
	if (iss.fail()) {
		cout<<filename<<"\n";
		cout<<"Can't extract time from file name "<<filename<<"!\n";
	    throw std::runtime_error{"failed to parse time string"};
	}
	return TDatime(tm.tm_year+1900,tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
}

vector<int> getFileLengthAndHeaders(TString filepath){
	ifstream file;
	file.open(filepath.Data());
	if (!file.is_open()){
		cout<<"cannot open "<<filepath<<"\n";
		throw std::runtime_error{"Cannot open file"};
	}
	char buff[256];
	file.getline(buff,256); //Header, names
	stringstream ss(buff);
	int Nvars = 0;
	int Apos = -1;
	int Bpos = -1;
	int Cpos = -1;
	int avgCpos = -1;
	int Nlines = 0;
	
	//Read headers
	while(ss.good()){
		string substr;
		getline(ss,substr,',');
		if (substr.find("average(C)") != string::npos){
			avgCpos = Nvars;
		} else if (substr.find("A") != string::npos){
			Apos = Nvars;
		} else if (substr.find("B") != string::npos){
			Bpos = Nvars;
		} else if (substr.find("C") != string::npos){
			Cpos = Nvars;
		}
		Nvars++;
	}
	file.getline(buff,256);
	file.getline(buff,256);
	
	//Count lines
	while (!file.eof()){
		file.getline(buff,256);
		if (file.eof()) break;
		Nlines++;
	}

	file.close();
	return {Nlines,Nvars,Apos,Bpos,Cpos,avgCpos};
}

vector<vector<double>> readFileTraces(TString filepath, int Nvars){
	ifstream file;
	file.open(filepath.Data());
	if (!file.is_open()){
		cout<<"cannot open "<<filepath<<"\n";
		throw std::runtime_error{"Cannot open file"};
	}
	char buff[256];
	file.getline(buff,256); //Header, names
	file.getline(buff,256); //Header, units
	file.getline(buff,256); //Header, blank row
	char comma;
	vector<vector<double>> vectors;
	vectors.resize(Nvars);
	while (!file.eof()){
		double var=0;
		for (int i=0; i<Nvars; i++){
			if (i>0) file>>comma;
			file>>var;
			if (file.eof()) break;
			vectors[i].push_back(var);
		}
	}
	file.close();
	return vectors;
}

void analyze_eddycurrents(TString folder, TString output_file, int Nfilesmax = -1){

	if (ROOT_VERSION_CODE < ROOT_VERSION(6,24,0)){
		cout<<"ERROR! You are using a ROOT version older than 6.24/06\n";
		cout<<"If you are on gm2ita, please execute:\n\n";
		cout<<"  source /opt/cernroot/bin/thisroot.sh\n\n";
		return;
	}

	double firstkick_max = 20.0;
	double kick_dt_guess = 10.0;
	double kick_trigger = -200.0;

	//Begin reading of files
	vector<TString> files = getListOfFiles(folder);
	int Nfiles = files.size();

	//Reading first file to get trace info
	vector<int> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
	int Nlines = headers[0];
	int Nvars = headers[1];

	float tstart = 0.;
	float tend = 100.;
	float dt = (tend-tstart)/Nlines;
	float t_before = -1.0;
	float t_after = 8.0;
	int Nkickbins = (t_after-t_before)/dt;
	float t_after8 = 22.0;
	int Nkickbins8 = (t_after8-t_before)/dt;

	TFile* fout = new TFile(output_file,"recreate");

	TProfile* g_fulltrace = new TProfile("trace","Trace;Time [ms];Voltage [mV]",Nlines,tstart,tend);
	TProfile** g_trace_kicks = new TProfile*[8];
	for (int i=0; i<8; i++){
		TString hname = Form("trace_kick%d",i+1);
		TString htitle = Form("Trace Kick %d;Time [ms];Voltage [mV]",i+1);
		g_trace_kicks[i] = new TProfile(hname,htitle,Nkickbins,t_before,t_after);
	}
	TProfile* g_trace_kick8long = new TProfile("trace_kick8long","Trace Kick 8;Time [ms];Voltage [mV]",Nkickbins8,t_before,t_after8);
	
	TProfile** g_trace_trend = new TProfile*[10];
	for (int i=0; i<10; i++){
		TString hname = Form("trace_trend_%d_%d",i*10,(i+1)*10);
		TString htitle = Form("Trace [%d-%d]%%;Time [ms];Voltage [mV]",i*10,(i+1)*10);
		g_trace_trend[i] = new TProfile(hname,htitle,Nlines,tstart,tend);
	}

	TProfile* g_trace_kickavg = (TProfile*)g_trace_kicks[0]->Clone("trace_avg");
	g_trace_kickavg->Reset();
	g_trace_kickavg->SetTitle("Kick average");

	TGraph* g_trend_A = new TGraph();
	TGraph* g_trend_B = new TGraph();
	TGraph* g_trend_ABdiff = new TGraph();
	TGraphErrors* g_trend_baseline = new TGraphErrors();
	TGraph* g_trend_stddev = new TGraph();
	TGraphErrors* g_trend_blumlein = new TGraphErrors();
	TGraph* g_trend_SNR = new TGraph();

	TGraph* g_correlation_AB_blumlein = new TGraph();
	TGraph* g_correlation_ABdiff_blumlein = new TGraph();
	TGraph* g_correlation_AB_SNR = new TGraph();
	TGraph* g_correlation_ABdiff_SNR = new TGraph();

	TGraph* g_blumlein_temp = new TGraph();
	TF1* f_blumlein = new TF1("f_blumlein","[0]+[2]*(x-[1])*(x-[1])",-0.5,0.0);
	double blumlein_fit_start = -0.4;
	double blumlein_fit_end = -0.2; 

	for (int fi=0; fi<Nfiles; fi++){

		if (Nfilesmax > 0 && fi >= Nfilesmax) break;

		TString fname = files[fi];
		if (fi%10==0) cout<<"Reading file \""<<fname<<"\" ("<<fi+1<<" of "<<Nfiles<<")\n";
		vector<pair<double,double>> trace;

		TDatime datetime = getFileTime(fname);

		TString filepath = Form("%s/%s",folder.Data(),fname.Data());

		//Reading file
		vector<vector<double>> traces = readFileTraces(filepath,Nvars);
		vector<double> trace_time = traces[0];
		vector<double> trace_A = traces[headers[2]];
		vector<double> trace_B = traces[headers[3]];
		vector<double> trace_C = traces[headers[4]];
		vector<double> trace_avgC = traces[headers[5]];
		
		if (trace_time.size() != Nlines){
			cout<<"Warning! File "<<fname<<" has "<<trace_time.size()<<" lines instead of the expected "<<Nlines<<"!\n";
		}

		//Find first kick
		double first_kick_guess = -1;
		for (int i=1; i<trace_time.size(); i++){
			if (trace_avgC[i-1] > kick_trigger && trace_avgC[i] < kick_trigger){
				first_kick_guess = trace_time[i];
				break;
			}
		}
		if (first_kick_guess < 0){
			cout<<"Warning! Could not find first kick in file "<<fname<<". Skipping...\n";
			continue;
		} else if (first_kick_guess > firstkick_max){
			cout<<"Warning! First kick found at "<<first_kick_guess<<" in file "<<fname<<". Skipping...\n";
			continue;
		}

		double A_avg = 0;
		double B_avg = 0;
		double avgC_baseline = 0;
		double avgC_stddev = 0;
		double avgC_squared = 0;
		double avgC_baselineerr = 0;
		int Nlines_baseline = 0;

		//Measure A, B, baseline, blumlein
		g_blumlein_temp->Set(0);
		for (int i=0; i<trace_time.size(); i++){

			A_avg += trace_A[i];
			B_avg += trace_B[i];

			if (trace_time[i] < first_kick_guess-0.5){
				avgC_baseline += trace_avgC[i];
				avgC_squared += trace_avgC[i]*trace_avgC[i];
				Nlines_baseline++;
			}

			if (trace_time[i] > first_kick_guess-0.5 && trace_time[i] < first_kick_guess){
				g_blumlein_temp->AddPoint(trace_time[i]-first_kick_guess,trace_avgC[i]);
			}
		}

		A_avg /= trace_time.size();
		B_avg /= trace_time.size();
		avgC_baseline /= Nlines_baseline;
		avgC_squared /= Nlines_baseline;
		avgC_stddev = sqrt(avgC_squared - avgC_baseline*avgC_baseline);
		avgC_baselineerr = avgC_stddev / sqrt(Nlines_baseline);

		f_blumlein->SetParameters(60.0,0.5*(blumlein_fit_start+blumlein_fit_end),-1000.0);
		TFitResultPtr fit_blumlein = g_blumlein_temp->Fit("f_blumlein","QS","",blumlein_fit_start,blumlein_fit_end);
		double peak = avgC_baseline;
		double peakerr = 0;
		if (fit_blumlein>=0){
			peak = fit_blumlein->Parameter(0);
			peakerr = fit_blumlein->ParError(0);
		}
		double blumlein = peak - avgC_baseline;
		double blumerr = sqrt(peakerr*peakerr + avgC_baselineerr*avgC_baselineerr);
		double AB = A_avg + B_avg;
		double ABdiff = B_avg - A_avg;
		double SNR = blumlein/avgC_stddev;

		g_trend_A->AddPoint(datetime.Convert(),A_avg);
		g_trend_B->AddPoint(datetime.Convert(),B_avg);
		g_trend_ABdiff->AddPoint(datetime.Convert(),ABdiff);
		g_trend_baseline->AddPoint(datetime.Convert(),avgC_baseline);
		g_trend_baseline->SetPointError(g_trend_baseline->GetN()-1,0,avgC_baselineerr);
		g_trend_stddev->AddPoint(datetime.Convert(),avgC_stddev);
		g_trend_blumlein->AddPoint(datetime.Convert(),blumlein);
		g_trend_blumlein->SetPointError(g_trend_blumlein->GetN()-1,0,blumerr);
		g_trend_SNR->AddPoint(datetime.Convert(),SNR);

		g_correlation_AB_blumlein->AddPoint(AB,blumlein);
		g_correlation_ABdiff_blumlein->AddPoint(ABdiff,blumlein);
		g_correlation_AB_SNR->AddPoint(AB,SNR);
		g_correlation_ABdiff_SNR->AddPoint(ABdiff,SNR);

		//Fill the full trace
		for (int i=0; i<trace_time.size(); i++){
			g_fulltrace->Fill(trace_time[i],trace_avgC[i] - avgC_baseline);

			int Nfi_10p = Nfiles / 10;
			int fi_p = (fi-fi%Nfi_10p)/Nfi_10p;
			if (fi_p<10){
				g_trace_trend[fi_p]->Fill(trace_time[i],trace_avgC[i] - avgC_baseline);
			}
		}

		for (int i=0; i<8; i++){
			double this_guess = first_kick_guess+i*kick_dt_guess-0.5;

			for (int j=0; j<trace_time.size(); j++){
				if (trace_time[j] > this_guess){
					if (trace_avgC[j-1] > kick_trigger){
						if (trace_avgC[j] < kick_trigger){
							double fraction = (kick_trigger-trace_avgC[j-1])/(trace_avgC[j]-trace_avgC[j-1]);
							double kick_time = trace_time[j-1] + dt*fraction;

							for (int k=0; k<trace_time.size(); k++){
								if (trace_time[k] >= kick_time + t_before && trace_time[k] < kick_time + t_after){
									g_trace_kicks[i]->Fill(trace_time[k]-kick_time,trace_avgC[k]);
								}
								if (i==7 && trace_time[k] >= kick_time + t_before && trace_time[k] < kick_time + t_after8){
									g_trace_kick8long->Fill(trace_time[k]-kick_time,trace_avgC[k]);
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
		g_trace_kickavg->Add(g_trace_kicks[i]);
	}


	g_trend_A->SetName("g_trend_A");
	g_trend_B->SetName("g_trend_B");
	g_trend_ABdiff->SetName("g_trend_ABdiff");
	g_trend_baseline->SetName("g_trend_baseline");
	g_trend_stddev->SetName("g_trend_stddev");
	g_trend_blumlein->SetName("g_trend_blumlein");
	g_trend_SNR->SetName("g_trend_SNR");
	g_correlation_AB_blumlein->SetName("g_correlation_AB_blumlein");
	g_correlation_ABdiff_blumlein->SetName("g_correlation_ABdiff_blumlein");
	g_correlation_AB_SNR->SetName("g_correlation_AB_SNR");
	g_correlation_ABdiff_SNR->SetName("g_correlation_ABdiff_SNR");

	g_trend_A->SetTitle("Trend of channel A;Time;A [V]");
	g_trend_B->SetTitle("Trend of channel B;Time;B [V]");
	g_trend_ABdiff->SetTitle("Trend of B-A;Time;B-A [V]");
	g_trend_baseline->SetTitle("Trend of avgC baseline;Time;Baseline [mV]");
	g_trend_stddev->SetTitle("Trend of avgC StdDev;Time;StdDev [mV]");
	g_trend_blumlein->SetTitle("Trend of blumlein peak;Time;Blumlein [mV]");
	g_trend_SNR->SetTitle("Trend of SNR;Time;Blumlein / StdDev");
	g_correlation_AB_blumlein->SetTitle("Correlation blumlein vs A+B;A+B [V];Blumlein [mV]");
	g_correlation_ABdiff_blumlein->SetTitle("Correlation blumlein vs B-A;B-A [V];Blumlein [mV]");
	g_correlation_AB_SNR->SetTitle("Correlation SNR vs A+B;A+B [V];Blumlein SNR");
	g_correlation_ABdiff_SNR->SetTitle("Correlation SNR vs B-A;B-A [V];Blumlein SNR");

	g_trend_A->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_A->GetXaxis()->SetTimeOffset(0);
	g_trend_A->GetXaxis()->SetTimeDisplay(1);
	g_trend_B->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_B->GetXaxis()->SetTimeOffset(0);
	g_trend_B->GetXaxis()->SetTimeDisplay(1);
	g_trend_ABdiff->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_ABdiff->GetXaxis()->SetTimeOffset(0);
	g_trend_ABdiff->GetXaxis()->SetTimeDisplay(1);
	g_trend_baseline->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_baseline->GetXaxis()->SetTimeOffset(0);
	g_trend_baseline->GetXaxis()->SetTimeDisplay(1);
	g_trend_stddev->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_stddev->GetXaxis()->SetTimeOffset(0);
	g_trend_stddev->GetXaxis()->SetTimeDisplay(1);
	g_trend_blumlein->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_blumlein->GetXaxis()->SetTimeOffset(0);
	g_trend_blumlein->GetXaxis()->SetTimeDisplay(1);
	g_trend_SNR->GetXaxis()->SetTimeFormat("%H:%M");
	g_trend_SNR->GetXaxis()->SetTimeOffset(0);
	g_trend_SNR->GetXaxis()->SetTimeDisplay(1);

	g_trend_A->Write();
	g_trend_B->Write();
	g_trend_ABdiff->Write();
	g_trend_baseline->Write();
	g_trend_stddev->Write();
	g_trend_blumlein->Write();
	g_trend_SNR->Write();

	g_correlation_AB_blumlein->Write();
	g_correlation_ABdiff_blumlein->Write();
	g_correlation_AB_SNR->Write();
	g_correlation_ABdiff_SNR->Write();

	fout->Write();
	fout->Close();
	cout<<"\nDone! Created output file "<<output_file<<"\n";
	return;
}

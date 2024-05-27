int nbins = 4692;
TH1D* EC_R0 = new TH1D("EC_R0","",nbins,0,700);
TH1D* EC_R0_cumulative = new TH1D("EC_R0_cumulative","",nbins,0,700);
TH1D* EC_R1 = new TH1D("EC_R1","",nbins,0,700);
TH1D* EC_R1_cumulative = new TH1D("EC_R1_cumulative","",nbins,0,700);

Double_t f_EC(Double_t* x, Double_t* par){
	double cum = EC_R0_cumulative->Interpolate(x[0]);
	return par[0]*exp(-x[0]/par[1])*(1+par[2]*cos(par[3]*(x[0]+cum)+par[4]));
}

double B0 = 1.45e7; // mG
double f_azimuth = 0.085;
double f_kickers = (53.1+53.0+55.0)/(3*55.0);
double mG_to_ppb = f_azimuth*f_kickers/B0;

double T_a_mus = 4.365233;
double tCyclotron = 0.14920;
double tau_gamma_mus_fake = 64.4;
std::pair<double, std::pair<double, double> > iN             = {1.0e+9, {1.e+5, 1.e+10}};
std::pair<double, std::pair<double, double> > itau           = {64.44, {63., 66.}};
std::pair<double, std::pair<double, double> > itau_ratio     = {tau_gamma_mus_fake, {tau_gamma_mus_fake, tau_gamma_mus_fake}};//Will get fixed value for Ratio
std::pair<double, std::pair<double, double> > iA             = {0.37, {0., 1.}};
std::pair<double, std::pair<double, double> > iR             = {1.43, {1.42, 1.44}};
std::pair<double, std::pair<double, double> > iphi           = {4.12, {0., 2*M_PI}};
std::pair<double, std::pair<double, double> > iA_CBO         = {0.05, {0., 1.}};
std::pair<double, std::pair<double, double> > iA_CBO_ratio   = {0.005, {0., 0.5}};//Smaller for ratio?
std::pair<double, std::pair<double, double> > iw_CBO         = {2.33, {2.1, 2.5}};
std::pair<double, std::pair<double, double> > ip_CBO         = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > itau_CBO       = {250., {50., 450.}};
std::pair<double, std::pair<double, double> > itau_CBO_ratio = {250., {150., 400.}};//Better like this for ratio?
std::pair<double, std::pair<double, double> > iA_VW          = {1.e-3, {0., 1.e-2}};
std::pair<double, std::pair<double, double> > ik_fy          = {1., {0.5, 1.5}};
std::pair<double, std::pair<double, double> > ip_VW          = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > ikLM           = {1.e-3, {-0.1, 0.1}};
std::pair<double, std::pair<double, double> > ikLM_ratio     = {1.e-3, {-0.5, 0.5}};//Better like this for ratio
std::pair<double, std::pair<double, double> > iA_CBOA        = {0.001, {0., 1.}};
std::pair<double, std::pair<double, double> > ip_CBOA        = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > iA_CBOP        = {0.001, {0., 1.}};
std::pair<double, std::pair<double, double> > ip_CBOP        = {0, {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > iA_2CBO        = {0.001, {0., 1.}};
std::pair<double, std::pair<double, double> > ip_2CBO        = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > itau_VW        = {45., {10., 150.}};
std::pair<double, std::pair<double, double> > iA_fy          = {0.0001, {0., 1.e-2}};
std::pair<double, std::pair<double, double> > ip_fy          = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > iw_VW          = {14.45, {13., 16.}};
std::pair<double, std::pair<double, double> > iA_cy          = {0.5, {0., 1.}};
std::pair<double, std::pair<double, double> > if_cy          = {6.70, {6., 7.4}};
std::pair<double, std::pair<double, double> > ip_cy          = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > itau_CBOVW     = {30., {2., 100.}};
std::pair<double, std::pair<double, double> > iA_CBOplusVW   = {0.0001, {0., 0.05}};
std::pair<double, std::pair<double, double> > ip_CBOplusVW   = {0., {-2*M_PI, 4*M_PI}};
std::pair<double, std::pair<double, double> > iA_CBOminusVW  = {0.001, {0., 0.05}};
std::pair<double, std::pair<double, double> > ip_CBOminusVW  = {0., {-2*M_PI, 4*M_PI}};

double waref = iR.first;

double Bk(double t){
    double cum = EC_R0->Integral(EC_R0->FindBin(0),EC_R0->FindBin(t));
    //EC_R0_cumulative->Interpolate(t);
    //if (t<30) cum = 0;
    return cum;
}

double Bkexp(double t){
	double A = -35 * mG_to_ppb;
	double tau = 47.4;
	return A*tau*(1-exp(-t/tau));
}

Double_t fitFunc12Par(Double_t *x, Double_t *par){
	Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    Double_t A = par[2];
    Double_t R = par[3];
    //limit phase to be between 0 and 2pi
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    Double_t ACBO = par[5];
    Double_t wCBO = par[6];
    //limit phase to be between 0 and 2pi
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
    Double_t tVW = tCBO * wCBO / wVW;
    Double_t w = R;
    Double_t CBOTerm = 1.0 + (exp(-t / tCBO) * ACBO * cos(wCBO * t - pCBO));
    Double_t VWTerm = 1.0 + exp(-t / tVW) * AVW * cos(wVW * t - pVW);
    return N * exp(-t / B) * (1 + A * cos( w * t - p)) * CBOTerm * VWTerm;
}

Double_t fitFunc12Par_EC(Double_t *x, Double_t *par){
    Double_t t = x[0];
    Double_t N = par[0];
    Double_t B = par[1];
    Double_t A = par[2];
    Double_t R = par[3];
    //limit phase to be between 0 and 2pi
    while (par[4] < 0)              par[4] += TMath::TwoPi();
    while (par[4] > TMath::TwoPi()) par[4] -= TMath::TwoPi();
    Double_t p = par[4];
    Double_t ACBO = par[5];
    Double_t wCBO = par[6];
    //limit phase to be between 0 and 2pi
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
    Double_t tVW = tCBO * wCBO / wVW;
    Double_t w = R;
    Double_t CBOTerm = 1.0 + (exp(-t / tCBO) * ACBO * cos(wCBO * t - pCBO));
    Double_t VWTerm = 1.0 + exp(-t / tVW) * AVW * cos(wVW * t - pVW);

    return N * exp(-t / B) * (1 + A * cos( w * (t+Bk(t)) - p)) * CBOTerm * VWTerm;
}

vector<double> pars_12 = {iN.first,itau.first,iA.first,iR.first,iphi.first,iA_CBO.first,iw_CBO.first,ip_CBO.first,itau_CBO.first,iA_VW.first,iw_VW.first,ip_VW.first};



void test_integral(){

	TFile* f_in_R0 = TFile::Open("../INFN_golden_R0_Bon.root");
	TGraph* g_EC_R0 = (TGraph*)f_in_R0->Get("g_R0_Bon");

	TFile* f_in_R1 = TFile::Open("../INFN_golden_R1_Bon.root");
	TGraph* g_EC_R1 = (TGraph*)f_in_R1->Get("g_R1_Bon");

	for (int bx=1; bx<=EC_R0->GetNbinsX(); bx++){
		double x = EC_R0->GetBinCenter(bx);
		double val = g_EC_R0->Eval(1e-3*x)*mG_to_ppb;
		//if (x<30) val = 0;
		//EC_R0->SetBinContent(bx,val);
		EC_R0->SetBinContent(bx,-35*mG_to_ppb*exp(-x/47.4));
	}
	for (int bx=1; bx<=EC_R1->GetNbinsX(); bx++){
		double x = EC_R1->GetBinCenter(bx);
		double val = g_EC_R1->Eval(1e-3*x)*mG_to_ppb;
		//if (x<30) val = 0;
		EC_R1->SetBinContent(bx,val);
	}
	EC_R0_cumulative = (TH1D*)EC_R0->GetCumulative();
	EC_R1_cumulative = (TH1D*)EC_R1->GetCumulative();

	TH1D* EC_analitic_cumulative = new TH1D("EC_analitic_cumulative","",nbins,0,700);
	for (int bx=1; bx<=EC_analitic_cumulative->GetNbinsX(); bx++){
		double x = EC_analitic_cumulative->GetBinCenter(bx);
		EC_analitic_cumulative->SetBinContent(bx,Bkexp(x));
	}

	new TCanvas();
	EC_R0_cumulative->Draw();
	EC_analitic_cumulative->SetLineColor(kRed);
	EC_analitic_cumulative->Draw("SAME");

	TF1* f_wiggle = new TF1("f_wiggle",fitFunc12Par,0,700,12);
	for (int i=0; i<12; i++){
		f_wiggle->SetParameter(i,pars_12[i]);
	}
    f_wiggle->SetNpx(2000);
    
	TF1* f_wiggle_EC = new TF1("f_wiggle_EC",fitFunc12Par_EC,0,700,12);
	for (int i=0; i<12; i++){
		f_wiggle_EC->SetParameter(i,pars_12[i]);
	}
    f_wiggle_EC->SetNpx(2000);

	TH1D* h1_wiggle = new TH1D("h1_wiggle","",nbins,0,700);
	TH1D* h1_wiggle_EC = new TH1D("h1_wiggle_EC","",nbins,0,700);
	TH1D* h1_wiggle_diff = new TH1D("h1_wiggle_diff","",nbins,0,700);

	for (int bx=1; bx<=h1_wiggle->GetNbinsX(); bx++){
		double t = h1_wiggle->GetBinCenter(bx);
		//if (t<30) continue;
		double N = f_wiggle->Eval(t);
		h1_wiggle->SetBinContent(bx,N);
		h1_wiggle->SetBinError(bx,sqrt(N));
	}
	for (int bx=1; bx<=h1_wiggle_EC->GetNbinsX(); bx++){
		double t = h1_wiggle_EC->GetBinCenter(bx);
		//if (t<30) continue;
		double N_EC = f_wiggle_EC->Eval(t);
		h1_wiggle_EC->SetBinContent(bx,N_EC);
		h1_wiggle_EC->SetBinError(bx,sqrt(N_EC));
	}
	for (int bx=1; bx<=h1_wiggle_diff->GetNbinsX(); bx++){
		double h1 = h1_wiggle->GetBinContent(bx);
		double h2 = h1_wiggle_EC->GetBinContent(bx);
		h1_wiggle_diff->SetBinContent(bx,h2-h1);
	}


	new TCanvas();
	EC_R0->Draw();
	EC_R1->Draw("SAME");

	new TCanvas();
	EC_R0_cumulative->Draw();
	EC_R1_cumulative->Draw("SAME");

	new TCanvas();
	h1_wiggle->Draw();
	for (int i=0; i<12; i++){
		f_wiggle->SetParameter(i,pars_12[i]);
	}
	h1_wiggle->Fit(f_wiggle,"","",30,650);
	double delta = 1e9*(f_wiggle->GetParameter(3)-waref)/waref;
	cout<<delta<<" ppb\n";

	new TCanvas();
	h1_wiggle_EC->Draw();
	for (int i=0; i<12; i++){
		f_wiggle->SetParameter(i,pars_12[i]);
	}
	h1_wiggle_EC->Fit(f_wiggle,"","",30,650);
	double delta_EC = 1e9*(f_wiggle->GetParameter(3)-waref)/waref;
	cout<<delta_EC<<" ppb\n";


	new TCanvas();
	h1_wiggle_diff->Draw("HIST");
}
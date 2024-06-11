{
	double binW = 0.1492;
	double xmin = 0;
	double xmax = 700;
	int nbins = (xmax-xmin)/binW;

	TH1D* h1_wiggle = new TH1D("h1_wiggle","",nbins,xmin,xmax);

	double w = 1.44;
	double A = 1;
	double phi0 = 0;

	double dw = (1-1e-7)*w;


	double T_end = xmax/2;

	TF1* f_wiggle = new TF1("f_wiggle","[0]*cos([1]*x+[2])",xmin,xmax);
	f_wiggle->SetNpx(nbins);
	f_wiggle->SetParNames("A","w","#phi");
	f_wiggle->SetParameters(A,w,phi0);
	
	TLine* l_A = new TLine(xmin,A,xmax,A);
	TLine* l_w = new TLine(xmin,w,xmax,w);
	TLine* l_phi0 = new TLine(xmin,phi0,xmax,phi0);

	//Scan for T_end
	TGraphErrors* g_tend_A = new TGraphErrors();
	TGraphErrors* g_tend_w = new TGraphErrors();
	TGraphErrors* g_tend_phi = new TGraphErrors();

	TH1D* h1_tend_1;
	TH1D* h1_tend_2;
	TH1D* h1_tend_3;

	new TCanvas();
	for (int i=0; i<nbins; i++){
		double T_end_i = i;
		h1_wiggle->Reset();
	
		double phit = phi0;
		double dx = h1_wiggle->GetBinWidth(1);
		for (int bx=1; bx<=h1_wiggle->GetNbinsX(); bx++){
			double x = h1_wiggle->GetBinCenter(bx);
			double wt = (bx<T_end_i?dw:w);
			if (bx==1) phit += wt*(dx*0.5);
			else phit += wt*dx;

			double val = A*cos(phit);
			h1_wiggle->SetBinContent(bx,val);
		}
		if (i==0){
			h1_tend_1 = (TH1D*)h1_wiggle->Clone("h1_tend_1");
		} else if (i==nbins/2){
			h1_tend_2 = (TH1D*)h1_wiggle->Clone("h1_tend_2");
		} else if (i==nbins-1){
			h1_tend_3 = (TH1D*)h1_wiggle->Clone("h1_tend_3");
		}

		f_wiggle->SetParameters(A,w,phi0);
		h1_wiggle->Draw("HIST");
		h1_wiggle->Fit(f_wiggle,"Q0+","",xmin,xmax);
		double fitted_A = f_wiggle->GetParameter(0);
		double fitted_A_err = f_wiggle->GetParError(0);
		double fitted_w = f_wiggle->GetParameter(1);
		double fitted_w_err = f_wiggle->GetParError(1);
		double fitted_phi = f_wiggle->GetParameter(2);
		double fitted_phi_err = f_wiggle->GetParError(2);
		double bias = (fitted_w-w)/w;

		g_tend_A->SetPoint(i,i,fitted_A);
		//g_tend_A->SetPointError(i,0,fitted_A_err);
		g_tend_w->SetPoint(i,i,fitted_w);
		//g_tend_w->SetPointError(i,0,fitted_w_err);
		g_tend_phi->SetPoint(i,i,fitted_phi);
		//g_tend_phi->SetPointError(i,0,fitted_phi_err);
	}

	g_tend_A->GetXaxis()->SetTitle("T_end");
	g_tend_w->GetXaxis()->SetTitle("T_end");
	g_tend_phi->GetXaxis()->SetTitle("T_end");
	g_tend_A->SetTitle("Fitted A");
	g_tend_w->SetTitle("Fitted w");
	g_tend_phi->SetTitle("Fitted #phi");


	TCanvas* can = new TCanvas("can","",1800,600);
	can->Divide(3,1);
	can->cd(1);
	//g_tend_A->GetYaxis()->SetRangeUser(0.8,1.2);
	g_tend_A->SetMarkerStyle(20);
	g_tend_A->Draw("APZ");
	l_A->Draw("SAME");
	can->cd(2);
	//g_tend_w->GetYaxis()->SetRangeUser(dw-(w-dw),w+(w-dw));
	g_tend_w->SetMarkerStyle(20);
	g_tend_w->Draw("APZ");
	l_w->Draw("SAME");
	can->cd(3);
	//g_tend_phi->GetYaxis()->SetRangeUser(-0.5,0.5);
	g_tend_phi->SetMarkerStyle(20);
	g_tend_phi->Draw("APZ");
	l_phi0->Draw("SAME");

	TCanvas* can_example = new TCanvas("can_example","",1800,600);
	can_example->Divide(3,1);
	can_example->cd(1);
	h1_tend_1->Draw();
	can_example->cd(2);
	h1_tend_2->Draw();
	can_example->cd(3);
	h1_tend_3->Draw();




	//Scan for dw
	TGraphErrors* g_dw_A = new TGraphErrors();
	TGraphErrors* g_dw_w = new TGraphErrors();
	TGraphErrors* g_dw_phi = new TGraphErrors();
	TH1D* h1_dw_1;
	TH1D* h1_dw_2;
	TH1D* h1_dw_3;
	new TCanvas();
	for (int i=0; i<1000; i++){
		double dw_i = (1+0.000001*(i-500))*w;
		h1_wiggle->Reset();
	
		double phit = phi0;
		double dx = h1_wiggle->GetBinWidth(1);
		for (int bx=1; bx<=h1_wiggle->GetNbinsX(); bx++){
			double x = h1_wiggle->GetBinCenter(bx);
			double wt = (bx<T_end?dw_i:w);
			if (bx==1) phit += wt*(dx*0.5);
			else phit += wt*dx;

			double val = A*cos(phit);
			h1_wiggle->SetBinContent(bx,val);
		}
		if (i==0){
			h1_dw_1 = (TH1D*)h1_wiggle->Clone("h1_dw_1");
		} else if (i==500){
			h1_dw_2 = (TH1D*)h1_wiggle->Clone("h1_dw_2");
		} else if (i==999){
			h1_dw_3 = (TH1D*)h1_wiggle->Clone("h1_dw_3");
		}

		f_wiggle->SetParameters(A,w,phi0);
		h1_wiggle->Draw("HIST");
		h1_wiggle->Fit(f_wiggle,"Q0+","",xmin,xmax);
		double fitted_A = abs(f_wiggle->GetParameter(0));
		double fitted_A_err = f_wiggle->GetParError(0);
		double fitted_w = f_wiggle->GetParameter(1);
		double fitted_w_err = f_wiggle->GetParError(1);
		double fitted_phi = fmod(f_wiggle->GetParameter(2),2*M_PI);
		double fitted_phi_err = f_wiggle->GetParError(2);
		double bias = (fitted_w-w)/w;

		g_dw_A->SetPoint(i,dw_i/w,fitted_A);
		//g_dw_A->SetPointError(i,0,fitted_A_err);
		g_dw_w->SetPoint(i,dw_i/w,fitted_w);
		//g_dw_w->SetPointError(i,0,fitted_w_err);
		g_dw_phi->SetPoint(i,dw_i/w,fitted_phi);
		//g_dw_phi->SetPointError(i,0,fitted_phi_err);
	}

	g_dw_A->GetXaxis()->SetTitle("dw");
	g_dw_w->GetXaxis()->SetTitle("dw");
	g_dw_phi->GetXaxis()->SetTitle("dw");
	g_dw_A->SetTitle("Fitted A");
	g_dw_w->SetTitle("Fitted w");
	g_dw_phi->SetTitle("Fitted #phi");


	TCanvas* can2 = new TCanvas("can2","",1800,600);
	can2->Divide(3,1);
	can2->cd(1);
	//g_dw_A->GetYaxis()->SetRangeUser(0.8,1.2);
	g_dw_A->SetMarkerStyle(20);
	g_dw_A->Draw("APZ");
	l_A->Draw("SAME");
	can2->cd(2);
	//g_dw_w->GetYaxis()->SetRangeUser(dw-(w-dw),w+(w-dw));
	g_dw_w->SetMarkerStyle(20);
	g_dw_w->Draw("APZ");
	l_w->Draw("SAME");
	can2->cd(3);
	//g_dw_phi->GetYaxis()->SetRangeUser(-0.5,0.5);
	g_dw_phi->SetMarkerStyle(20);
	g_dw_phi->Draw("APZ");
	l_phi0->Draw("SAME");

	TCanvas* can2_example = new TCanvas("can2_example","",1800,600);
	can2_example->Divide(3,1);
	can2_example->cd(1);
	h1_dw_1->Draw();
	can2_example->cd(2);
	h1_dw_2->Draw();
	can2_example->cd(3);
	h1_dw_3->Draw();

}
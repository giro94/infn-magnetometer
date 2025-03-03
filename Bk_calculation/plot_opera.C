{

    TGraph *g_opera = new TGraph();
    int i=0;

    g_opera->SetPoint(i, -45.1, 0); i++;
    g_opera->SetPoint(i, -42.73969179419197, 65.14270675169547); i++;
    g_opera->SetPoint(i, -41.63201281902776, 66.39766170665047); i++;
    g_opera->SetPoint(i, -40.5058432114496, 67.50480817896548); i++;
    g_opera->SetPoint(i, -39.218792231360284, 68.5423625873064); i++;
    g_opera->SetPoint(i, -37.529537819993046, 69.55461079056585); i++;
    g_opera->SetPoint(i, -35.77995289393413, 70.5636957181901); i++;
    g_opera->SetPoint(i, -33.970037453183517, 71.43675979350137); i++;
    g_opera->SetPoint(i, -32.039460983049537, 71.6929851199514); i++;
    g_opera->SetPoint(i, -30.108884512915557, 71.76890373519586); i++;
    g_opera->SetPoint(i, -28.178308042781572, 71.44624962040692); i++;
    g_opera->SetPoint(i, -26.247731572647592, 71.13308533252354); i++;
    g_opera->SetPoint(i, -24.43781613189698, 70.46879744913453); i++;
    g_opera->SetPoint(i, -22.56757017645469, 69.73175422613625); i++;
    g_opera->SetPoint(i, -20.817985250395767, 68.74164895232312); i++;
    g_opera->SetPoint(i, -19.00806980964516, 67.82113574248405); i++;
    g_opera->SetPoint(i, -17.19815436889455, 66.93858184026723); i++;
    g_opera->SetPoint(i, -15.388238928143947, 66.07500759186152); i++;
    g_opera->SetPoint(i, -13.578323487393336, 65.23041299726692); i++;
    g_opera->SetPoint(i, -11.768408046642724, 64.42377771029456); i++;
    g_opera->SetPoint(i, -09.898162091200433, 63.649829604885774); i++;
    g_opera->SetPoint(i, -08.088246650449822, 62.89907885413502); i++;
    g_opera->SetPoint(i, -06.218000695007531, 62.41193440631642); i++;
    g_opera->SetPoint(i, -04.2874242248735506, 62.165198906771934); i++;
    g_opera->SetPoint(i, -02.3568477547395617, 62.0608108108108); i++;
    g_opera->SetPoint(i, -00.4262712846055816, 62.01336167628301); i++;
    g_opera->SetPoint(i, 01.5646357002200872, 62.01336167628301); i++;
    g_opera->SetPoint(i, 03.555542685045756, 62.00387184937746); i++;
    g_opera->SetPoint(i, 05.365458125796367, 62.35499544488307); i++;
    g_opera->SetPoint(i, 07.296034595930347, 62.772547828727596); i++;
    g_opera->SetPoint(i, 09.226611066064327, 63.3039781354388); i++;
    g_opera->SetPoint(i, 11.096857021506619, 63.81010223706852); i++;
    g_opera->SetPoint(i, 12.90677246225723, 64.72112562000201); i++;
    g_opera->SetPoint(i, 14.656357388316152, 65.60051624658365); i++;
    g_opera->SetPoint(i, 16.526603343758444, 66.43351216384923); i++;
    g_opera->SetPoint(i, 18.336518784509055, 67.31922934170125); i++;
    g_opera->SetPoint(i, 20.146434225259657, 68.1290279043088); i++;
    g_opera->SetPoint(i, 21.89601915131858, 69.16869116307316); i++;
    g_opera->SetPoint(i, 23.70593459206919, 70.07971454600667); i++;
    g_opera->SetPoint(i, 25.515850032819802, 70.89583965988459); i++;
    g_opera->SetPoint(i, 27.325765473570414, 71.56961737017915); i++;
    g_opera->SetPoint(i, 29.256341943704385, 72.00614940783478); i++;
    g_opera->SetPoint(i, 31.186918413838374, 72.07257819617368); i++;
    g_opera->SetPoint(i, 33.177825398664034, 71.82373384620574); i++;
    g_opera->SetPoint(i, 34.927410324722956, 71.00022775584571); i++;
    g_opera->SetPoint(i, 36.737325765473567, 70.05124506528999); i++;
    g_opera->SetPoint(i, 38.4265801768408, 68.86501670209535); i++;
    g_opera->SetPoint(i, 40.05550407351635, 67.56385599082226); i++;
    g_opera->SetPoint(i, 41.62409745550022, 66.15514390795288); i++;
    g_opera->SetPoint(i, 43.13258096010324, 64.51311150636988); i++;
    g_opera->SetPoint(i, 45.1, 0); i++;
    	
    
    new TCanvas();
    g_opera->SetMarkerStyle(20);
    g_opera->Draw("APL");


    //BIOT SAVART

   	TH1D* h1_x0 = new TH1D("h1_x0","By (x=0);y [mm];By",100,-50,50);
	TH1D* h1_xdR = new TH1D("h1_xdR","By (x=dR);y [mm];By",100,-50,50);
	TH1D* h1_y0 = new TH1D("h1_y0","By (y=0);x [mm];By",100,-50,50);


	TH2D* h2_xy = new TH2D("h2_xy","Field from many wires;x [mm];y [mm];By",100,-50,50,100,-50,50);
	TGraph* g_wires = new TGraph();

	double D = 91.4;
	double R = D/2.;
	double maxTheta = 31. * (3.14159/180.);
	double minr = 1.0;
	double Bmax = 20;
	double dR = 17.50;


	double w_crystal = 4;
	double L_crystal = 32;
	TGraph* g_crystal1 = new TGraph();
	TGraph* g_crystal2 = new TGraph();
	g_crystal1->SetPoint(0,-w_crystal/2,-L_crystal/2);
	g_crystal1->SetPoint(1,-w_crystal/2,L_crystal/2);
	g_crystal1->SetPoint(2,w_crystal/2,L_crystal/2);
	g_crystal1->SetPoint(3,w_crystal/2,-L_crystal/2);
	g_crystal1->SetPoint(4,-w_crystal/2,-L_crystal/2);

	g_crystal2->SetPoint(0,dR-w_crystal/2,-L_crystal/2);
	g_crystal2->SetPoint(1,dR-w_crystal/2,L_crystal/2);
	g_crystal2->SetPoint(2,dR+w_crystal/2,L_crystal/2);
	g_crystal2->SetPoint(3,dR+w_crystal/2,-L_crystal/2);
	g_crystal2->SetPoint(4,dR-w_crystal/2,-L_crystal/2);

	bool first = true;
	for (double x=-R; x<R; x+=1){
		for (double y=-R; y<R; y+=1){
			double B=0;
			for (double a=-maxTheta; a<maxTheta; a+=0.001){

				double w1x = -R*cos(a);
				double w1y = R*sin(a);
				double w2x = R*cos(a);
				double w2y = R*sin(a);
				if (first){
					g_wires->AddPoint(w1x,w1y);
					g_wires->AddPoint(w2x,w2y);
				}

				double r1 = sqrt((x-w1x)*(x-w1x) + (y-w1y)*(y-w1y));
				double r2 = sqrt((x-w2x)*(x-w2x) + (y-w2y)*(y-w2y));
				double By1 = 1/r1;
				double By2 = 1/r2;

				if (r1<minr) By1 = 0;
				if (r2<minr) By2 = 0;

				B += 0.1*By1;
				B += 0.1*By2;

			}
			first = false;

			if (x*x+y*y > R*R) continue;
		
			h2_xy->Fill(x,y,B);

			if (abs(x)<0.5){
				h1_x0->Fill(y,B);
			}
			if (abs(x-dR)<0.5){
				h1_xdR->Fill(y,B);
			}
			if (abs(y)<0.5){
				h1_y0->Fill(x,B);
			}

		}
	}

	cout<<"Using "<<g_wires->GetN()<<" wires\n";

	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);
	h2_xy->SetContour(64);

	new TCanvas("","",800,800);
	g_wires->SetMarkerStyle(20);
	g_wires->SetMarkerColor(kBlack);
	g_crystal1->SetLineWidth(2);
	g_crystal2->SetLineWidth(2);
	g_crystal1->SetMarkerColor(kRed);
	g_crystal2->SetMarkerColor(kRed);
	h2_xy->Draw("COLZ");
	g_wires->Draw("P");
	g_crystal1->Draw("L");
	g_crystal2->Draw("L");

	new TCanvas();
	h1_y0->SetLineWidth(2);
	h1_y0->Draw("HIST");

	new TCanvas();
	h1_x0->GetYaxis()->SetRangeUser(3.5,5.5);
	h1_xdR->GetYaxis()->SetRangeUser(3.5,5.5);
	h1_x0->SetLineWidth(2);
	h1_xdR->SetLineWidth(2);
	h1_x0->SetLineColor(kBlue);
	h1_xdR->SetLineColor(kRed);
	h1_x0->Draw("HIST");
	h1_xdR->Draw("HIST SAME");
	gPad->BuildLegend();


	double B0 = h1_y0->Interpolate(0);
	double BdR = h1_y0->Interpolate(dR);
	cout<<"B at x=0: "<<B0<<"\n";
	cout<<"B at x=dR: "<<BdR<<"\n";
	cout<<"Ratio : "<<BdR/B0<<"\n";


	TFitResultPtr res0 = h1_x0->Fit("pol0","S","",-L_crystal/2,L_crystal/2);
	TFitResultPtr resdR = h1_xdR->Fit("pol0","S","",-L_crystal/2,L_crystal/2);

	double By_R0 = res0->Parameter(0);
	double By_R1 = resdR->Parameter(0);
	cout<<"By in R0 crystal (fit) : "<<By_R0<<"\n";
	cout<<"By in R1 crystal (fit) : "<<By_R1<<"\n";
	cout<<"Ratio : "<<By_R1/By_R0<<"\n";

	double By_R0_int = h1_x0->Integral(h1_x0->FindBin(-L_crystal/2),h1_x0->FindBin(L_crystal/2));
	double By_R1_int = h1_xdR->Integral(h1_xdR->FindBin(-L_crystal/2),h1_xdR->FindBin(L_crystal/2));
	cout<<"By in R0 crystal (integral) : "<<By_R0_int<<"\n";
	cout<<"By in R1 crystal (integral) : "<<By_R1_int<<"\n";
	cout<<"Ratio : "<<By_R1_int/By_R0_int<<"\n";





	h1_y0->Scale(g_opera->Eval(0)/h1_y0->Interpolate(0));

	new TCanvas();
	h1_y0->Draw("HIST");
	g_opera->Draw("PL");


	TH1D* h1_diff = (TH1D*)h1_y0->Clone("h1_diff");
	for (int i=1; i<=h1_diff->GetNbinsX(); i++){
		double x = h1_diff->GetBinCenter(i);
		double y = h1_diff->GetBinContent(i);
		
		h1_diff->SetBinContent(i,y-g_opera->Eval(x));
	}

	new TCanvas();
	h1_diff->Draw("HIST");


	TF1* f0_4 = new TF1("f0_4","[0]+[1]*x*x+[2]*x^4+[3]*x^6",-45,45);

	f0_4->SetParameters(0,-0.07529,4.11e-5,1e-6);
	f0_4->FixParameter(0,0);
	//f0_4->FixParameter(1,0);
	f0_4->FixParameter(3,0);

	h1_diff->Fit(f0_4,"","",-45,45);
	f0_4->Draw("SAME");






}
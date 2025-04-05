{

	TFile* fin = TFile::Open("UMass_model_0p1.root");

	TH2D* h2 = (TH2D*)fin->Get("h2");
	TH2D* h1_x0 = (TH2D*)fin->Get("h1_x0");
	TH2D* h1_x17p5 = (TH2D*)fin->Get("h1_x17p5");
	TH2D* h1_y0 = (TH2D*)fin->Get("h1_y0");

	double B00 = h2->Interpolate(0,0);
	double B10 = h2->Interpolate(17.5,0);

	double r0 = 0.0;
	double r1 = 17.5;
	double dx = 2.0;
	double dy = 2.0;
	double ymin = -16;
	double ymax = 16;

	double yminUMASS = -14.5;
	double ymaxUMASS = 14.5;

	double r0bin = h2->GetXaxis()->FindBin(r0);
	double r1bin = h2->GetXaxis()->FindBin(r1);
	double yminbin = h2->GetYaxis()->FindBin(ymin);
	double ymaxbin = h2->GetYaxis()->FindBin(ymax);
	double R0_avg = h2->Integral(r0bin,r0bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);
	double R1_avg = h2->Integral(r1bin,r1bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);


	double yminbinUMASS = h2->GetYaxis()->FindBin(yminUMASS);
	double ymaxbinUMASS = h2->GetYaxis()->FindBin(ymaxUMASS);
	double R0_avgUMASS = h2->Integral(r0bin,r0bin,yminbinUMASS,ymaxbinUMASS) / (1+ymaxbinUMASS-yminbinUMASS);
	double R1_avgUMASS = h2->Integral(r1bin,r1bin,yminbinUMASS,ymaxbinUMASS) / (1+ymaxbinUMASS-yminbinUMASS);

	cout<<"Field at (0,0): "<<B00<<"\n";
	cout<<"Field at (17.5,0): "<<B10<<"\n";
	cout<<"Ratio R1/R0: "<<B10/B00<<"\n";

	cout<<"Measured at R0: "<<R0_avg<<" (norm factor: "<<R0_avg/B00<<")\n";
	cout<<"Measured at R1: "<<R1_avg<<" (norm factor: "<<R1_avg/B10<<")\n";
	cout<<"Ratio R1/R0: "<<R1_avg/R0_avg<<"\n";

	cout<<"UMass Measured at R0: "<<R0_avgUMASS<<" (norm factor: "<<R0_avgUMASS/B00<<")\n";
	cout<<"UMass Measured at R1: "<<R1_avgUMASS<<" (norm factor: "<<R1_avgUMASS/B10<<")\n";
	cout<<"Ratio R1/R0: "<<R1_avgUMASS/R0_avgUMASS<<"\n";

	//x +- 2 mm
	cout<<"\nx +- 2 mm\n";

	double B00p2 = h2->Interpolate(0+dx,0);
	double B10p2 = h2->Interpolate(17.5+dx,0);
	double r0p2bin = h2->GetXaxis()->FindBin(r0+dx);
	double r1p2bin = h2->GetXaxis()->FindBin(r1+dx);
	double R0p2_avg = h2->Integral(r0p2bin,r0p2bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);
	double R1p2_avg = h2->Integral(r1p2bin,r1p2bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);
	cout<<"Measured at R0+2: "<<R0p2_avg<<" (norm factor: "<<R0p2_avg/B00p2<<")\n";
	cout<<"Measured at R1+2: "<<R1p2_avg<<" (norm factor: "<<R1p2_avg/B10p2<<")\n";
	cout<<"Ratio R1/R0: "<<R1p2_avg/R0p2_avg<<"\n";
	cout<<"Ratio factor R0p2/R0: "<<R0p2_avg/R0_avg<<"\n";

	double B00m2 = h2->Interpolate(0-dx,0);
	double B10m2 = h2->Interpolate(17.5-dx,0);
	double r0m2bin = h2->GetXaxis()->FindBin(r0-dx);
	double r1m2bin = h2->GetXaxis()->FindBin(r1-dx);
	double R0m2_avg = h2->Integral(r0m2bin,r0m2bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);
	double R1m2_avg = h2->Integral(r1m2bin,r1m2bin,yminbin,ymaxbin) / (1+ymaxbin-yminbin);
	cout<<"Measured at R0-2: "<<R0m2_avg<<" (norm factor: "<<R0m2_avg/B00m2<<")\n";
	cout<<"Measured at R1-2: "<<R1m2_avg<<" (norm factor: "<<R1m2_avg/B10m2<<")\n";
	cout<<"Ratio R1/R0: "<<R1m2_avg/R0m2_avg<<"\n";
	cout<<"Ratio factor R0m2/R0: "<<R0m2_avg/R0_avg<<"\n";



	//y +- 2 mm
	cout<<"\ny +- 2 mm\n";

	double B00yp2 = h2->Interpolate(0,0+dy);
	double B10yp2 = h2->Interpolate(17.5,0+dy);
	double yminp2bin = h2->GetYaxis()->FindBin(ymin+dy);
	double ymaxp2bin = h2->GetYaxis()->FindBin(ymax+dy);
	double R0yp2_avg = h2->Integral(r0bin,r0bin,yminp2bin,ymaxp2bin) / (1+ymaxp2bin-yminp2bin);
	double R1yp2_avg = h2->Integral(r1bin,r1bin,yminp2bin,ymaxp2bin) / (1+ymaxp2bin-yminp2bin);
	cout<<"Measured at R0 y+2: "<<R0yp2_avg<<" (norm factor: "<<R0yp2_avg/B00yp2<<")\n";
	cout<<"Measured at R1 y+2: "<<R1yp2_avg<<" (norm factor: "<<R1yp2_avg/B10yp2<<")\n";
	cout<<"Ratio R1/R0: "<<R1yp2_avg/R0yp2_avg<<"\n";
	cout<<"Ratio factor R0yp2/R0: "<<R0yp2_avg/R0_avg<<"\n";

	double B00ym2 = h2->Interpolate(0,0-dy);
	double B10ym2 = h2->Interpolate(17.5,0-dy);
	double yminm2bin = h2->GetYaxis()->FindBin(ymin-dy);
	double ymaxm2bin = h2->GetYaxis()->FindBin(ymax-dy);
	double R0ym2_avg = h2->Integral(r0bin,r0bin,yminm2bin,ymaxm2bin) / (1+ymaxm2bin-yminm2bin);
	double R1ym2_avg = h2->Integral(r1bin,r1bin,yminm2bin,ymaxm2bin) / (1+ymaxm2bin-yminm2bin);
	cout<<"Measured at R0 y-2: "<<R0ym2_avg<<" (norm factor: "<<R0ym2_avg/B00ym2<<")\n";
	cout<<"Measured at R1 y-2: "<<R1ym2_avg<<" (norm factor: "<<R1ym2_avg/B10ym2<<")\n";
	cout<<"Ratio R1/R0: "<<R1ym2_avg/R0ym2_avg<<"\n";
	cout<<"Ratio factor R0ym2/R0: "<<R0ym2_avg/R0_avg<<"\n";


	//y +- 5 mm
	cout<<"\ny +- 5 mm\n";

	double B00yp5 = h2->Interpolate(0,0+5);
	double B10yp5 = h2->Interpolate(17.5,0+5);
	double yminp5bin = h2->GetYaxis()->FindBin(ymin+5);
	double ymaxp5bin = h2->GetYaxis()->FindBin(ymax+5);
	double R0yp5_avg = h2->Integral(r0bin,r0bin,yminp5bin,ymaxp5bin) / (1+ymaxp5bin-yminp5bin);
	double R1yp5_avg = h2->Integral(r1bin,r1bin,yminp5bin,ymaxp5bin) / (1+ymaxp5bin-yminp5bin);
	cout<<"Measured at R0 y+2: "<<R0yp5_avg<<" (norm factor: "<<R0yp5_avg/B00yp5<<")\n";
	cout<<"Measured at R1 y+2: "<<R1yp5_avg<<" (norm factor: "<<R1yp5_avg/B10yp5<<")\n";
	cout<<"Ratio R1/R0: "<<R1yp5_avg/R0yp5_avg<<"\n";
	cout<<"Ratio factor R0yp5/R0: "<<R0yp5_avg/R0_avg<<"\n";

	double B00ym5 = h2->Interpolate(0,0-5);
	double B10ym5 = h2->Interpolate(17.5,0-5);
	double yminm5bin = h2->GetYaxis()->FindBin(ymin-5);
	double ymaxm5bin = h2->GetYaxis()->FindBin(ymax-5);
	double R0ym5_avg = h2->Integral(r0bin,r0bin,yminm5bin,ymaxm5bin) / (1+ymaxm5bin-yminm5bin);
	double R1ym5_avg = h2->Integral(r1bin,r1bin,yminm5bin,ymaxm5bin) / (1+ymaxm5bin-yminm5bin);
	cout<<"Measured at R0 y-2: "<<R0ym5_avg<<" (norm factor: "<<R0ym5_avg/B00ym5<<")\n";
	cout<<"Measured at R1 y-2: "<<R1ym5_avg<<" (norm factor: "<<R1ym5_avg/B10ym5<<")\n";
	cout<<"Ratio R1/R0: "<<R1ym5_avg/R0ym5_avg<<"\n";
	cout<<"Ratio factor R0ym5/R0: "<<R0ym5_avg/R0_avg<<"\n";

}
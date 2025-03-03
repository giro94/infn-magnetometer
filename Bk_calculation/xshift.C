{

	double R0 = 0;
	double R1 = 17.5;
	double B0 = -16.2;
	double B1 = -35.5;

	double dx = -2;

	double p1 = (-B1+B0)/((R0-dx)*(R0-dx) - (R1-dx)*(R1-dx));
	double p0 = B0 - p1*(R0-dx)*(R0-dx);

	cout<<"dx: "<<dx<<" p0: "<<p0<<" p1: "<<p1<<"\n";


}
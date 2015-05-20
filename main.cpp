#include "main.h"
#include "revision.h"

class Simulation
{
    // Working variables and output
    double frontCryst, oldCryst; // position of crystallization front as fractional index [-]
    double frontConv, oldConv; // position of convecting front as fractional index [-]
    double Qcmb; // CMB heat flow [W]
    double Qtop; // heat flow through the TBL [W]
    double nQsec, oQsec; // secular heat flow in the BMO [W]
    double nQlat, oQlat; // latent heat release [W]
    double *Q; // current heat distribution [J/m3]
    double *Q2; // previous heat distribution [J/m3]
    double *K, *P, *S;
    double Tcore, liquidEnrichment;

    // getters
    double getS(int); // solidus temperature
    double getT(int); // temperature
    double getP(int); // solid fraction
    double getK(int); // diffusivity
    double getC(int); // specific heat
    double getRadio(double); // heat production rate in w/kg
    
    // boundary tracking
    double func(double);
    double bissect(double);
    double getConvectiveFront(double);
    double getCrystallizationFront(void);

    // fluxes tracking
    void updateHeat(void);
    void getHeatContent(void);
    
    // numerical scheme
    double gradTF(int);
    double gradTB(int);
    double ThermalDiffusion(int);
    
    // logic
    void initialize(void);
    void iterate(double);
    void printInfo(void);

    
public:
    void run();


};

// default
double D = 4000e3;
double densityDrop = 1000;
double liquidusDrop = 1000;

int XR = (D+delta)/dx;

double lastTime = 0;
bool finished = false;
int lastBdy = 100000;
int bdyMoved = 1;
bool endProg = false;
bool convective = false;

string fname1, fname2, fname3;

/*--------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 --------------------------------------------------------------------------------*/

double Simulation::getRadio(double time)
// input time is in s after 4.5 Ga ago
// returns production rate in W/kg
{
    double t = 4.5e9 - time/3600./24/365;
    double H = 0.9928*CU0*HU238*exp(t*log(2.)/tU238);
    H += 0.0071*CU0*HU235*exp(t*log(2.)/tU235);
    H += CTh0*HTh232*exp(t*log(2.)/tTh232);
    H += 1.19e-4*CK0*HK40*exp(t*log(2.)/tK40);
	double prefactor = 1;
	if (conserveHeat == true)
	{
		prefactor = (pow(RC+D,3)-pow(RC,3))/(pow(RC+frontCryst*dx,3)-pow(RC,3));
	}
    return prefactor*H;
}

double Simulation::getS(int x)
// returns solidus temperature at cell x
{
    double T = TL0 - liquidusDrop*(D-x*dx)/D;
	double gradTL = 1000; // K between c=0 and c=1 (assumption for test purposes, check later)
	if (convective) gradTL = 2000;
	double dc = 0.0;

	if(x < int(CMB/dx)) return 0.0;

	if(frontConv == int(CMB/dx) and x == XR-1) // this should be done only once per timestep
	{
		dc = 3*pow(RC+frontCryst*dx,2)*eta*(oldCryst-frontCryst)*dx;
		dc /= (pow(RC+frontCryst*dx,3)-pow(RC+frontConv*dx,3));
		liquidEnrichment += dc;
		if(liquidEnrichment > 0.5) liquidEnrichment = 0.5;
	}

	// the TBL starts solid, liquidus is irrelevant
    if(x*dx >= D) return 9999; 

	// the convection zone has the same liquidus as its lowermost point
    if(x > frontConv and convective) return TL0 - liquidEnrichment*gradTL;
    if(x > frontConv) return TL0 - liquidusDrop*(D-frontConv*dx)/D - liquidEnrichment*gradTL;

    return T - P[x]*liquidusDrop*dx/D;
}


double Simulation::getT(int x)
// returns temperature at cell x
{
    if (x < 0) return Tcore;
    if (x > XR-1) return TMantle;
    
    double T = S[x];
    double Qv = Q[x] - T*C;
    
    if (Qv < 0) T = Q[x]/C; // below Ts
    if (Qv > Qlatent) T = S[x] + (Qv - Qlatent)/C; // above Ts
    
    return T;
}

double Simulation::getP(int x)
// phi is solid fraction
{
    if (getT(x) > S[x]) return 0;
    if (getT(x) < S[x]) return 1;
    
    if (x < 0) return 0;
    if (x > XR-1) return 1;
    
    return 1-(Q[x]-C*S[x])/Qlatent;
}

double Simulation::func(double r)
// frontCryst as a function of frontConv
// this is derived analytically (cf. manuscript)
// and is inverted to get frontConv knowing frontCryst.
{
    double s = densityDrop/drho/D; // [wt.%/m]
    double B = eta/s;
    double rtop = RC + D;
    double u, g;
    
    if(r == rtop) return r;
    
    u = exp(-(rtop-r)/B);
    
    g = 6*pow(B,3)*(1-u);
    g += 6*pow(B,2)*(r-rtop*u);
    g += 3*B*(pow(r,2)-pow(rtop,2)*u);
    g += pow(r,3);
    g = pow(g, 1./3);
    
    return g;
}

double Simulation::bissect(double target)
{
    // bissection algorithm to invert "func"
    int i=0, imax = 100;
    double tol = 5e-1;
    double a=RC, b=RC+frontCryst*dx, c;
    double ca, cb, cc;
    
    assert(target < RC+D);
    assert(func(RC+D) == RC+D);
    
    if (func(a) > target) return RC;
    
    while (i <= imax){
        c = (a+b)/2;
        ca = func(a)-target;
        cb = func(b)-target;
        cc = func(c)-target;
        if ( cc == 0 or (b-a)/2 < tol) {return c;}
        if ( signbit(cc) == signbit(ca) ) { a = c; } else { b = c; }
        i++;
    }
    assert(i < imax); // otherwise we haven't converged
    return 0.0;
}

double Simulation::getConvectiveFront(double x)
// x is the floating index corresponding to frontCryst
// here we need to find r such that g(r) = x*dx + RC
{
    if (finished == true or convective) return int(CMB/dx);
    double target = RC+x*dx;
    double r = target;
    
    if (target == RC+D) return D/dx;
    assert(target < RC+D);
    
    r = bissect(target);
    double rconv = (r-RC)/dx;
    
    if (rconv < int(CMB/dx)){ finished = true; return int(CMB/dx); }
    return rconv;
}


double Simulation::getCrystallizationFront(void)
// get crystallization front position (as fraction)
// from solid fraction vector
// Nich: Force forward motion
{
    int i = XR;
    for (int j=XR; j>=0; j--)
    {
        if (P[j] > 0) { i = j; }
    }
    
    if (i > lastBdy) i = lastBdy;
    
    if (i != lastBdy)
    {
		bdyMoved = 1;
		lastBdy = i;
	}
	
	if (i == int(CMB/dx) and P[i] == 1.) { bdyMoved = 1; endProg = true; }
    if (i == XR) { return i;}
    return i + (1-P[i]);
}

void Simulation::updateHeat(void)
{
    Qcmb += 4*PI*pow(RC,2)*C*K[0]*(Tcore-getT(0))/dx*dt; // J

    int x = frontCryst;
    
    double ke = 2*K[x+1]*K[x]/(K[x+1]+K[x]);
    double qe = ke*gradTF(x);
    
    Qtop += 4*PI*pow(RC+frontCryst*dx,2)*qe*C*dt; // J
}

void Simulation::getHeatContent(void)
{
    oQlat= nQlat;
    nQlat = 0;
    oQsec = nQsec;
    nQsec = 0;
    
    for(int i=0; i<XR; i++)
    {
        double stability = 0.5*(gradTF(i)+gradTB(i))*alpha*rho - densityDrop/D*dx;
        if (stability > 0 and i < frontCryst)
        {
            ofstream logfile;
            logfile.open (fname1.c_str());
            logfile << i << " should convect..." << endl;
            logfile.close();
        }
        
        nQlat += 4*PI*pow(RC+i*dx,2)*dx*Qlatent*(1-P[i]); // J
        nQsec += 4*PI*pow(RC+i*dx,2)*dx*(getT(i)-S[i])*C*(1-P[i]); // J
    }
}

double Simulation::getK(int x)
// get diffusivity as a function of position
{
    // convecting mantle
    if (x >= frontCryst + TBL) {return KTC;}
    
    // within TBL the diffusivity should be KTM
    if ((x >= frontCryst) && (x < frontCryst + TBL)) {return KTM;}

    // everywhere else should be KTL
    return KTL;
}

double Simulation::getC(int x)
// get specific heat as a function of position
// only one value so far; kept for possible future improvement
{
    double phi = P[x];
    return phi*C + (1-phi)*C;
}

// div K * grad(T)
double Simulation::gradTF(int x) {return (getT(x)-getT(x+1))/dx;}
double Simulation::gradTB(int x) {return (getT(x-1)-getT(x))/dx;}

double Simulation::ThermalDiffusion(int x)
// now in spherical coordinates
{
    double r = RC + x*dx;
    
    double kw, ke;
    double qw, qe, fw, fe;
    double rw, re;
    
    rw = r - 0.5*dx;
    re = r + 0.5*dx;
    
    kw = 2*K[x-1]*K[x]/(K[x-1]+K[x]);
    ke = 2*K[x+1]*K[x]/(K[x+1]+K[x]);
    
    double gTB = gradTB(x);
    double gTF = gradTF(x);
    qw = kw*gTB;
    qe = ke*gTF;

    fw = 0;
    fe = 0;
    
    double Tx = getT(x);
    
    // within convective zone of the BMO, should be KTS
    if ((x <= floor(frontCryst) - 1) && (x > ceil(frontConv)) + 1)
    {
        
        double cp = 1000;
        double TaB = 0.5*(getT(x-1)+Tx);
        double TaF = 0.5*(getT(x+1)+Tx);
        
        double dTB = max(0., gTB - alpha*g*TaB/cp);
        double dTF = max(0., gTF - alpha*g*TaF/cp);
        
        double mu = 1e16/rho;
        double LB = min(frontCryst-x+0.5,x-frontConv-0.5)*dx;
        double LF = min(frontCryst-x-0.5,x-frontConv+0.5)*dx;
        if (LF < 0) LF = 0;
        if (LB < 0) LB = 0;
        
        double KcB = min(1e-3, alpha*g*pow(LB,4)/(18*mu)*dTB);
        double KcF = min(1e-3, alpha*g*pow(LF,4)/(18*mu)*dTF);
        
        fw = KcB*dTB;
        fe = KcF*dTF;
    }

    if (x == XR-1) {fe = 0; qe = K[x]*(Tx-TMantle)/dx;}
    if (x == 0) {fw = 0; qw = K[x]*(Tcore-Tx)/dx;}
    
    // equivalent to the commented one
    return ((qw+fw)*pow(rw,2) - (qe+fe)*pow(re,2))/dx/pow(r,2);
//    return (qw*pow(rw,2) - qe*pow(re,2))/dx/pow(r,2) + (fw*pow(rw,2) - fe*pow(re,2))/dx/pow(r,2);
}

void Simulation::iterate(double time)
{
    double *buf;
    double dQ, S0;
    
    oldCryst = frontCryst;
    frontCryst = getCrystallizationFront();
    oldConv = frontConv;
    frontConv = getConvectiveFront(frontCryst);

	if (time == 0) { oldCryst = frontCryst; oldConv = frontConv; }
    
    assert(! isnan(frontCryst));
    assert(! isnan(frontConv));
    
    for (int x=0;x<XR;x++)
    {
        dQ = dt*ThermalDiffusion(x)*getC(x) + dt*rho*getRadio(time);
        assert(! isnan(dQ));

        Q2[x] = Q[x] + dQ;
        K[x] = getK(x);
        P[x] = getP(x);
		S0 = S[x];
        S[x] = getS(x);
		assert(! isnan(S[x]));
		// artificial change of liquidus with time needs to be taken into
		// account to keep energy balance consistent
        nQsec += 4*PI*pow(RC+x*dx,2)*dx*(S0-S[x])*C*(1-P[x]); // J
    }
    
    Tcore = Tcore - 3/RC*KTL*dt*(Tcore-getT(0))/dx; // W/m3 * s * W/m/K
    
    buf = Q; Q = Q2; Q2 = buf;
    
    return;
}

void Simulation::initialize()
{
    Q=(double*)malloc(XR*sizeof(double));
    Q2=(double*)malloc(XR*sizeof(double));
    K=(double*)malloc(XR*sizeof(double));
    P=(double*)malloc(XR*sizeof(double));
    S=(double*)malloc(XR*sizeof(double));
    
    Tcore = TL0 + dT0*D/1000.; // 1 K / km of initial gradient
	liquidEnrichment = 0;

    frontCryst = XR; //getCrystallizationFront();
    frontConv = XR; //getConvectiveFront(frontCryst);

    // isotherm at Tcore, gradient in TBL
    for (int x=0;x<XR;x++)
    {
        S[x] = getS(x);
        
        if(x*dx < D)
        {
            Q[x] = C*S[x] + Qlatent + (Tcore-S[x]-dT0*x)*C; // 1K/km of initial gradient
        } else {
            Q[x] = C*(Tcore - (x*dx-D)*(Tcore-TMantle)/delta);
        }
        
        P[x] = getP(x);
        K[x] = getK(x);

    }
    
    getHeatContent();
    oQlat = nQlat;
    oQsec = nQsec;
    
    oldConv = frontConv;
    oldCryst = frontCryst;
    
    printInfo();
}

void Simulation::run()
{
    initialize();
    
    double time = 0;
    double oldtime = 0;
    int lastOut = 0;
    
    while (time < tmax and frontCryst >= frontConv and endProg != true)
    {
        iterate(time);
        updateHeat();
        time += dt;
     
        if (bdyMoved or time == dt)
        {
            getHeatContent();
            double massBMO = 4*PI*rho*(pow(RC+dx*frontCryst,3) - pow(RC,3))/3;
            
            if (time == 0) time = 1;

            FILE *g = fopen(fname3.c_str(), "a");
            fprintf(g, "%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
                    time/Ma,
                    (RC+dx*frontCryst)/1e3,
                    (RC+dx*frontConv)/1e3,
                    Qcmb/(time-lastTime)/1e12,
                    Qtop/(time-lastTime)/1e12,
                    (oQsec-nQsec)/(time-oldtime)/1e12,
                    (oQlat-nQlat)/(time-oldtime)/1e12,
                    getRadio(time)*massBMO/1e12);
            fclose(g);
            lastTime = time;
            oldtime = time;
            
            Qcmb = 0;
            Qtop = 0;
            bdyMoved = 0;
            assert(Qtop >= 0);
        }
        
        if (time >= snap*(1+lastOut))
        {
            ofstream logfile;
            logfile.open (fname1.c_str(), ios::app);
            logfile << "t = " << setprecision(4) << time/Ma << " Ma" << endl;
            logfile.close();
            
            FILE *f = fopen(fname2.c_str(), "a");
            for (int x=0;x<XR;x++)
            fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
                    (RC+x*dx)/1e3,
                    time/Ma,
                    Q[x],
                    getT(x),
                    S[x],
                    P[x],
                    K[x]);
            
            fprintf(f,"\n");
            fclose(f);
            lastOut++;
        }
    }
}

void Simulation::printInfo(void)
{
    ofstream logfile;
    logfile.open (fname1.c_str());
    
    logfile << "revision = " << revision << endl;
    logfile << "dx = " << dx << " m" << endl;
    logfile << "dt = " << dt/Ma << " Ma" << endl;
    logfile << "N-factor = " << snap/dt << " (n-timesteps per snapshots)" << endl;
    logfile << "tmax = " << tmax/Ma << " Ma" << endl;
    logfile << "nsteps = " << ceil(tmax/dt) << endl;
    logfile << "ncells = " << XR << endl;
    logfile << endl;
    
    logfile << "KTL = " << KTL << endl;
    logfile << "KTM = " << KTM << endl;
    logfile << "KTC = " << KTC << endl;
    logfile << "rho = " << rho << endl;
    logfile << "Qlatent = " << Qlatent << endl;
    logfile << "C = " << C << endl;
    logfile << "alpha = " << alpha << endl;
    logfile << "g = " << g << endl;
    logfile << "drho = " << drho << endl;
    logfile << "densityDrop = " << densityDrop << endl;
    logfile << "liquidusDrop = " << liquidusDrop << endl;
    logfile << "eta = " << eta << endl;
    logfile << "dT0 = " << dT0 << endl;
    logfile << "TMantle = " << TMantle << endl;
    logfile << "TL0 = " << TL0 << endl;
    
    logfile.close();
}

int main(int argc, char **argv)
{
    if (argc != 5) {cout << "Wrong number of arguments!" << endl; exit(1);}

	D = atof(argv[1])*1000;
	liquidusDrop = atof(argv[2]);
	densityDrop = atof(argv[3]);
	XR = (D+delta)/dx;
	if (atoi(argv[4]) == 1) convective = true;

	ostringstream str; 
	str << "D" << D/1000 << "L" << liquidusDrop << "R" << densityDrop;
	string prepend = str.str();

    fname1 = prepend + string(".log");
    fname2 = prepend + string("-record.txt");
    fname3 = prepend + string("-timeseries.txt");
    
    Simulation s;
    s.run();
}

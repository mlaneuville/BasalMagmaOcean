#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

#include <cmath>
#include <cassert>

using namespace std;

class Simulation
{
    // useful declarations
    double PI = 3.14159265359;
    double years = 365*24*3600;
    double Ma = 1e6*years;

    // Model geometry
    double RC = 3480e3; // Core Radius [m]
    double VC = 4*PI*pow(RC,3)/3; // Core Volume [m3]
    double D = 200e3; // BMO thickness [m]
    double dx = 1.0e3; // Grid resolution [m]
    double delta = 100e3; // Mantle thermal boundary layer (TBL) [m]
    double CMB = 25e3; // CMB thermal boundary layer thickness [m]
    
    int TBL = delta/dx; // number of cells in TBL [-]
    int XR = (D+delta)/dx; // total number of cells [-]

    // Thermo-chemical properties
    double KTL = 1e-5; // Liquid diffusivity [m2/s]
    double KTM = 2e-6; // Mantle diffusivity through TBL [m2/s]
    double KTC = 1e-3; // Convecting mantle effective diffusivity [m2/s]
    double rho = 5000; // average BMO density [kg/m3]
    double Qlatent = rho*5000*300; // latent heat of crystallization [J/m3]
    double C = rho*1000; // specific heat [J/K/m3]
    double alpha = 2e-5; // thermal expansivity [/K]
    double g = 11.0; // CMB gravity [m/s2]
    
    /* @TODO: input should be in kg/m3/km and implications should be
              dealt with in func() and getS(). */
    double drho = 2.5e6; // "density" gradient in the BMO [m/wt.%]
                         // drho = D/dc, where dc is the compo difference
                         // between bottom and top of the BMO
    double eta = 0.088;  // Fe partitioning upon crystallization [wt.%]
    
    // Numerical variables
    double tmax = 150*Ma; // max simulation time [Ma]
    double snap = 0.1*Ma; // snapshot frequency [Ma]
    double dt = 0.25*dx*dx/KTC; // simulation timestep [s]

    // initial conditinos
    double Tcore = 4500; // initial core temperature [K]
    double TMantle = 2500; // initial mantle temperature [K]

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

    // getters
    double getS(int); // solidus temperature
    double getT(int); // temperature
    double getP(int); // solid fraction
    double getK(int); // diffusivity
    double getC(int); // specific heat

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
    void iterate(void);

public:
    void run();

};

double lastTime = 0;
bool finished = false;
int lastBdy = 100000;
int bdyMoved = 1;

/*--------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------
 --------------------------------------------------------------------------------*/

double Simulation::getS(int x)
// returns solidus temperature at cell x
{
    double TA = 3500; // solidus temperature at the base of the BMO [K]
    double TB = 4500; // solidus temperature at the top of the BMO [K]
    double T = TB - (TB-TA)*(D-(x+1)*dx)/D;
    if(x*dx>=D) return 5000; // the TBL starts solid, liquidus is irrelevant
    return T - P[x]*(TB-TA)*dx/D;
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
// this is derived analytically (C. manuscript)
// and is inverted to get frontConv knowing frontCryst.
{
    double B = eta*drho;
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
    if (finished == true) return int(CMB/dx);
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
    if ((x <= floor(frontCryst)) && (x > ceil(frontConv)))
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
        
        double KcB = alpha*g*pow(LB,4)/(18*mu);
        double KcF = alpha*g*pow(LF,4)/(18*mu);
        
        fw = KcB*pow(dTB,2);
        fe = KcF*pow(dTF,2);
    }

    if (x == XR-1) {fe = 0; qe = K[x]*(Tx-TMantle)/dx;}
    if (x == 0) {fw = 0; qw = K[x]*(Tcore-Tx)/dx;}
    
    // equivalent to the commented one
    return ((qw+fw)*pow(rw,2) - (qe+fe)*pow(re,2))/dx/pow(r,2);
//    return (qw*pow(rw,2) - qe*pow(re,2))/dx/pow(r,2) + (fw*pow(rw,2) - fe*pow(re,2))/dx/pow(r,2);
}

void Simulation::iterate()
{
    double *buf;
    double dQ;
    
    oldCryst = frontCryst;
    frontCryst = getCrystallizationFront();
    oldConv = frontConv;
    frontConv = getConvectiveFront(frontCryst);
    
    assert(! isnan(frontCryst));
    assert(! isnan(frontConv));
    
    for (int x=0;x<XR;x++)
    {
        dQ = dt*ThermalDiffusion(x)*getC(x);
        assert(! isnan(dQ));

        Q2[x] = Q[x] + dQ;
        K[x] = getK(x);
        P[x] = getP(x);
        S[x] = getS(x);
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
    
    frontCryst = XR; //getCrystallizationFront();
    frontConv = XR; //getConvectiveFront(frontCryst);
    
    // isotherm at Tcore, gradient in TBL
    for (int x=0;x<XR;x++)
    {
        S[x] = getS(x);
        
        if(x*dx < D)
        {
            Q[x] = C*S[x] + Qlatent + (Tcore-S[x])*C;
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
    
    cout << "dx = " << dx << " m" << endl;
    cout << "dt = " << dt/Ma << " Ma" << endl;
    cout << "N-factor = " << snap/dt << " (n-timesteps per snapshots)" << endl;
    cout << "tmax = " << tmax/Ma << " Ma" << endl;
    cout << "nsteps = " << ceil(tmax/dt) << endl;
    cout << "ncells = " << XR << endl;
    cout << endl;

}

void Simulation::run()
{
    initialize();
    
    double time = 0;
    double oldtime = 0;
    int lastOut = 0;
    
    while (time < tmax)
    {
        iterate();
        updateHeat();
        time += dt;
     
        if (bdyMoved)
        {
            getHeatContent();
            
            if (time == 0) time = 1;
            assert(Qtop >= 0);

            FILE *g = fopen("timeseries.txt", "a");
            fprintf(g, "%.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
                    time/Ma,
                    (RC+dx*frontCryst)/1e3,
                    (RC+dx*frontConv)/1e3,
                    Qcmb/(time-lastTime)/1e12,
                    Qtop/(time-lastTime)/1e12,
                    (oQsec-nQsec)/(time-oldtime)/1e12,
                    (oQlat-nQlat)/(time-oldtime)/1e12);
            fclose(g);
            lastTime = time;
            oldtime = time;
            
            Qcmb = 0;
            Qtop = 0;
            bdyMoved = 0;
        }
        
        if (time >= snap*(1+lastOut))
        {
            cout << "t = " << setprecision(4) << time/Ma << " Ma" << endl;
            FILE *f = fopen("record.txt", "a");
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

int main(int argc, char **argv)
{
    Simulation s;
    s.run();
}

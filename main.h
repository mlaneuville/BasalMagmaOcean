#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

// useful declarations
static const double PI = 3.15;
static const double years = 365*24*3600;
static const double Ma = 1e6*years;

// Model geometry
static const double RC = 3480e3; // Core Radius [m]
static const double VC = 4*PI*pow(RC,3)/3; // Core Volume [m3]
static const double D = 200e3; // BMO thickness [m]
static const double dx = 1.0e3; // Grid resolution [m]
static const double delta = 100e3; // Mantle thermal boundary layer (TBL) [m]
static const double CMB = 25e3; // CMB thermal boundary layer thickness [m]

static const int TBL = delta/dx; // number of cells in TBL [-]
static const int XR = (D+delta)/dx; // total number of cells [-]

// Thermo-chemical properties
static const double KTL = 1e-5; // Liquid diffusivity [m2/s]
static const double KTM = 2e-6; // Mantle diffusivity through TBL [m2/s]
static const double KTC = 1e-3; // Convecting mantle effective diffusivity [m2/s]
static const double rho = 5000; // average BMO density [kg/m3]
static const double Qlatent = rho*5000*300; // latent heat of crystallization [J/m3]
static const double C = rho*1000; // specific heat [J/K/m3]
static const double alpha = 2e-5; // thermal expansivity [/K]
static const double g = 11.0; // CMB gravity [m/s2]
    
/* @TODO: input should be in kg/m3/km and implications should be
          dealt with in func() and getS(). */
static const double drho = 2.5e6; // "density" gradient in the BMO [m/wt.%]
                     // drho = D/dc, where dc is the compo difference
                     // between bottom and top of the BMO
static const double eta = 0.088;  // Fe partitioning upon crystallization [wt.%]

// Numerical variables
static const double tmax = 1000*Ma; // max simulation time [Ma]
static const double snap = 0.1*Ma; // snapshot frequency [Ma]

// initial conditinos
static const double TMantle = 2500; // initial mantle temperature [K]

static const double dt = 0.25*dx*dx/KTC; // simulation timestep [s]



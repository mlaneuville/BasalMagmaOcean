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
static const double PI = 3.14159265359;
static const double years = 365*24*3600;
static const double Ma = 1e6*years;

// Model geometry
static const double RC = 3480e3; // Core Radius [m]
static const double VC = 4*PI*pow(RC,3)/3; // Core Volume [m3]
static const double dx = 1.0e3; // Grid resolution [m]
static const double delta = 100e3; // Mantle thermal boundary layer (TBL) [m]
static const double CMB = 25e3; // CMB thermal boundary layer thickness [m]

static const int TBL = delta/dx; // number of cells in TBL [-]

// Thermo-chemical properties
static const double KTL = 1e-5; // Liquid diffusivity [m2/s]
static const double KTM = 2e-6; // Mantle diffusivity through TBL [m2/s]
static const double KTC = 1e-3; // Convecting mantle effective diffusivity [m2/s]
static const double rho = 5000; // average BMO density [kg/m3]
static const double Qlatent = rho*5000*300; // latent heat of crystallization [J/m3]
static const double C = rho*1000; // specific heat [J/K/m3]
static const double alpha = 2e-5; // thermal expansivity [/K]
static const double g = 11.0; // CMB gravity [m/s2]
    
// Compositional gradient definition
static const double drho = 6e3; // density gradient [kg/m3/wt.%] - no need to change this
static const double eta = 0.088;  // Fe partitioning upon crystallization [wt.%]

// Numerical variables
static const double tmax = 4500*Ma; // max simulation time [Ma]
static const double snap = 1*Ma; // snapshot frequency [Ma]
static const double dt = 0.25*dx*dx/KTC; // simulation timestep [s]

// initial conditions
static const bool conserveHeat = true;
static const double dT0 = 1.; // initial temperature gradient in BMO [K/km]
static const double TL0 = 4500.; // liquidus at the top of the layer [K]

// Specific heat production rate [W/kg]
static const double HU238 = 9.46e-5;
static const double HU235 = 5.69e-4;
static const double HTh232 = 2.64e-5;
static const double HK40 = 2.92e-5;

// Half life of the isotope [a]
static const double tU238 = 4.47e9;
static const double tU235 = 7.04e8;
static const double tTh232 = 1.40e10;
static const double tK40 = 1.25e9;

// Initial concentrations
static const double CU0 = 20.3e-9;
static const double CTh0 = 79.5e-9;
static const double CK0 = 240e-6;


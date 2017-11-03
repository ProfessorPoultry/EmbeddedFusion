#include "FusionCalculator.h"
#include "FusionCrossSectionSolver.h"
#include "Maxwellian.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

static const double k 		 	= 1.38064852e-23;
static const double m_proton 	= 1.6726219e-27;
static const double eVtoJ 		= 1.60218e-19;
void main() {
	FusionCalculator FuseCalc;
	FusionCrossSectionSolver FuseSolver;
	Maxwellian Maxwellian;
	
	int  	Emin 			= 1000;
	int 	Emax			= 3000;
	int 	Erange 			= 1000;
	double 	Ue 				= 800; 
	double  n1				= 1e16;
	double	n2 				= 6.76e22; //density of palladium surface 
	// double  n2				= 1e14; //density of tokamak

	//Tritium
	double 	m1 				= 2 * m_proton;
	double 	m2 				= 3 * m_proton;
	double 	Z1 				= 1;
	double 	Z2 				= 1;
	int 	type			= 2; //Tritium
	FuseCalc.EnergyRange(m1, m2, Z1, Z2, Emin, Emax, Erange, Ue, type);


	// double  n2				= 1;
	// double 	Ue 				= 800; 
	// double 	rate 			= FuseCalc.ReactionRate(n1, n2, m1, m2, Z1, Z2, E, Ue, SvarArray);
	// double 	energyOutTrit 	= rate/2*(1010000 + 3020000)*eVtoJ;
	// double 	energyOutHe 	= rate/2*(820000 +  2450000)*eVtoJ;
	// double 	energyOut 		= energyOutTrit + energyOutHe;
	// cout << "energyOutt << endl; 
	// double energyIn = E * n2 * eVtoJ;
	// cout << "energyIn " << energyIn << endl;
}
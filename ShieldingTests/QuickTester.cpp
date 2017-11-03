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
	Maxwellian Maxwell;
	//deuterium
	double 	m1 				= 3 * m_proton;
	double 	m2 				= 2 * m_proton;
	double 	Z1 				= 1;
	double 	Z2 				= 1;
	int 	type 			= 2;
	int 	E 				= 3000;
	double	n1 				= 1e16; 
	double  n2				= 1.66e16;//density of palladium surface 
	double 	Ue 				= 800; 
	FuseSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, type);
	double vr = Maxwell.ConvertKineticEnergyToSpeed(E*eVtoJ, FuseSolver.mr);
	cout << "vr = " << vr << endl;
	double rate = n1*n2*FuseSolver.P*vr;
	cout << "rate = " << rate << endl;
	double 	energyOutTrit 	= rate*(3.5e6 + 14.1e6);
	double 	energyOut 		= energyOutTrit;

	cout << "energyOut " << energyOut << endl; 
	double energyIn = E * n1;
	cout << "energyIn " << energyIn << endl;
}
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
	FusionCrossSectionSolver FuseSolver;
	Maxwellian Maxwell;
	//deuterium
	double 	m1 				= 3 * m_proton;
	double 	m2 				= 2 * m_proton;
	double 	Z1 				= 1;
	double 	Z2 				= 1;
	int 	type 			= 2;

	int 	E 				= 2000;
	double  n1				= 1e15; //ion flux
	double	n2 				= 6.76e22; //density of palladium surface 
	double 	energyOutTrit 	= (3.5e6 + 14.1e6);

	double 	vr 				= 0;
	double 	rate			= 0;
	ofstream crossSections ("UeRange/crossSections.txt");
	ofstream Ues ("UeRange/Energies.txt");
	ofstream InputEnergy ("UeRange/InputEnergy.txt");
	ofstream outPutEnergy ("UeRange/outPutEnergy.txt");
	ofstream Rates ("UeRange/Rates.txt");
	for (int Ue = 27; Ue < 2000; Ue++)
	{
		Ues 			<< Ue << endl;
		FuseSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, type);
		crossSections 	<< FuseSolver.P << endl; 
		vr = Maxwell.ConvertKineticEnergyToSpeed(E*eVtoJ, FuseSolver.mr);
		rate = n1*n2*FuseSolver.P*vr;
		Rates << rate << endl;
		outPutEnergy << energyOutTrit * rate << endl;
		InputEnergy << n1*E << endl;

	}
}
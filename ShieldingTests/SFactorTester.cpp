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
	//deuterium
	double 	m1 				= 2 * m_proton;
	double 	m2 				= m1;
	double 	Z1 				= 1;
	double 	Z2 				= 1;
	std::vector<int> 		SvarArray;
	SvarArray.push_back(0);
	SvarArray.push_back(1);
	ofstream Sfactors ("SFactor/Sfactors.txt");
	ofstream Energies ("SFactor/Energies.txt");
	for (int E = 5000; E < 20000; E=E+1)
	{
		FuseSolver.calculateReducedMass(m1, m2);
		FuseSolver.calculateReducedEnergy(m1, E);
		FuseSolver.calculateSFactor(SvarArray);
		Sfactors << FuseSolver.S << endl;
		Energies << FuseSolver.Er << endl;
	}
}
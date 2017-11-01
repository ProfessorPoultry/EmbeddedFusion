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
	//deuterium
	double 	m1 				= 2 * m_proton;
	double 	m2 				= m1;
	double 	Z1 				= 1;
	double 	Z2 				= 1;
	std::vector<int> 		SvarArray;
	SvarArray.push_back(0);
	SvarArray.push_back(1);
	int 	E 				= 3000;
	double	n1 				= 1; 
	double  n2				= 1.66e15;//density of palladium surface 
	double 	Ue 				= 600; 
	double 	rate 			= FuseCalc.ReactionRate(n1, n2, m1, m2, Z1, Z2, E, Ue, SvarArray);
	double 	energyOutTrit 	= rate/2*(1010000 + 3020000);
	double 	energyOutHe 	= rate/2*(820000 +  2450000);
	double 	energyOut 		= energyOutTrit + energyOutHe;
	cout << "energyOut " << energyOut << endl; 
	double energyIn = E * n1;
	cout << "energyIn " << energyIn << endl;
}
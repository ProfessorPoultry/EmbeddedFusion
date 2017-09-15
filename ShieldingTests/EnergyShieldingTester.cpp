#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "FusionCrossSectionSolver.h"

using namespace std;

static const double k 		 	= 1.38064852e-23;
static const double m_proton 	= 1.6726219e-27;
static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;

void SingleEnergy (double m1, double m2, double Z1, double Z2) {
	double E 						= 1000 * eVtoJ;
	double Ue						= 27 *eVtoJ;
	double fusionProb 				= FCSSolver.returnFusionCross(E, m1, m2, Z1, Z2);
	double shieldedFusionProb		= FCSSolver.returnFusionCross(E+Ue, m1, m2, Z1, Z2);
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Unshielded: " 			<< fusionProb << endl;
	cout << "Shielded: " 			<< shieldedFusionProb << endl;
	cout << "E/Ue: " 				<< E/Ue << endl;
	cout << "EnhancementFactor: " 	<< enhancementFactor << endl;
}

void EnergyRange (double m1, double m2, double Z1, double Z2) {
	int Emin						= 1000; 	// in eV
	int Emax						= 3000; //in eV
	double Ue						= 400 *eVtoJ;
	int Esize						= Emax-Emin;
	double fusionProb 				[Esize];
	double shieldedFusionProb		[Esize];
	double enhancementFactors		[Esize];
	double Ej;
	int Eindex;
	ofstream UnshieldedProbs ("UnshieldedProbs.txt");
	ofstream ShieldedProbs ("ShieldedProbs.txt");
	ofstream Energies ("Energies.txt");
	ofstream EnhancementFactors ("EnhancementFactors.txt");
	for (int E = Emin; E <= Emax; ++E)
	{
		Ej 							= E * eVtoJ;
		Eindex = E - Emin;
		fusionProb[Eindex]			= FCSSolver.returnFusionCross(Ej, m1, m2, Z1, Z2);
		shieldedFusionProb[Eindex] 	= FCSSolver.returnFusionCross(Ej+Ue, m1, m2, Z1, Z2);
		enhancementFactors[Eindex]	= shieldedFusionProb[Eindex]/fusionProb[Eindex];
		UnshieldedProbs 			<< fusionProb[Eindex] << endl;
		ShieldedProbs				<< shieldedFusionProb[Eindex] << endl;
		EnhancementFactors			<< enhancementFactors[Eindex] << endl;
		Energies					<< Ej << endl;
	}

}

void main () {	
	double m1					= 2 * m_proton;
	double m2 					= m1;
	int Z1 = 1;
	int	Z2 = 1;
	EnergyRange(m1, m2, Z1, Z2);
}


#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "FusionCrossSectionSolver.h"

using namespace std;

static const double k 		 	= 1.38064852e-23;
static const double m_proton 	= 1.6726219e-27;
static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;

void SingleEnergy (double m1, double m2, double Z1, double Z2, double E, double Ue) {
	double Ej 						= E * eVtoJ;
	double Uej						= Ue *eVtoJ;
	double fusionProb 				= FCSSolver.returnFusionCross(Ej, m1, m2, Z1, Z2);
	double shieldedFusionProb		= FCSSolver.returnFusionCross(Ej+Uej, m1, m2, Z1, Z2);
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Unshielded: " 			<< fusionProb << endl;
	cout << "Shielded: " 			<< shieldedFusionProb << endl;
	cout << "E/Ue: " 				<< Ej/Uej << endl;
	cout << "EnhancementFactor: " 	<< enhancementFactor << endl;
}

void EnergyRange (double m1, double m2, double Z1, double Z2, int Emin, int Emax, double Ue) {
	double Uej						= Ue *eVtoJ;
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
		shieldedFusionProb[Eindex] 	= FCSSolver.returnFusionCross(Ej+Uej, m1, m2, Z1, Z2);
		enhancementFactors[Eindex]	= shieldedFusionProb[Eindex]/fusionProb[Eindex];
		UnshieldedProbs 			<< fusionProb[Eindex] << endl;
		ShieldedProbs				<< shieldedFusionProb[Eindex] << endl;
		EnhancementFactors			<< enhancementFactors[Eindex] << endl;
		Energies					<< Ej << endl;
	}
	cout << "Done, see text files in ShieldingTests";
}

void main () {	
	double m1;
	double m2;
	int Z1;
	int	Z2;

	string preset = "";
	int presetType;
	int type;

	while(type != 1 && type != 2 && type != 3) {
		cout 	<< "Would you like to use: " 
				<< endl << "1. a single particle energy" 
				<< endl << "2. a range" 
				<< endl << "3. Default Preset?" << endl;
		cin >> type;
		if (type != 1 && type != 2 && type != 3)
		{
			cout << "error, invalid response" << endl;
		}
	}
	if (type == 3)
	{
		cout << "D-D Fusion, E=1:3keV, Ue = 27eV" << endl;
		EnergyRange(2, 2, 1, 1, 1000,3000, 25);
		return 0;
	}

	cout << "would you like to use a particle preset? y/n" << endl;
	cin >> preset;
	if (preset == "y")
	{
		while (presetType != 1 && presetType != 2) {
			cout << "Presets:" << endl << "1. D-D fusion" << endl << "2. Cancel" << endl;
			cin >> presetType;
			if (presetType == 1)
			{
				m1 = 2 * m_proton;
				m2 = m1;
				Z1 = 1;
				Z2 = 1;
			}
			else if(presetType == 2){				
				break;
			}
			else{
				cout << "error, invalid response" << endl;
			}		
		}
	}
	else{
		cout << "Please enter m1, m2 in amu and Z1, Z2 as integers" << endl;
		cout << "m1" << endl;
		cin >> m1;
		cout << "m2" << endl;
		cin >> m2;
		cout << "Z1" << endl;
		cin >> Z1;
		cout << "Z2" << endl;
		cin >> Z2;
	}	
	if (type == 1)
	{
		int E;
		double Ue;
		cout << "please enter E in eV" << endl;
		cin >> E;
		cout << "please enter Ue in eV" << endl;
		cin >> Ue;
		SingleEnergy(m1, m2, Z1, Z2, E, Ue);
	}
	else if (type == 2)
	{
		int Emin;
		int Emax;
		double Ue;
		cout << "please enter Emin in eV" << endl;
		cin >> Emin;
		cout << "please enter Emax in eV" << endl;
		cin >> Emax;
		cout << "please enter Ue in eV" << endl;
		cin >> Ue;
		EnergyRange(m1, m2, Z1, Z2, Emin, Emax, Ue);
	}
	else {
		cout << "ERROR" <<endl;
		return 0;
	}	
}


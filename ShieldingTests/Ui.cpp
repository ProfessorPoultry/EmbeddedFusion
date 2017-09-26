
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

main () {	
	FusionCalculator FuseCalc;
	FusionCrossSectionSolver FuseSolver;

	double m1;
	double m2;
	
	int Z1;
	int	Z2;

	string preset = "";
	int presetType;
	int type;
	std::vector<int> SvarArray;

	while(!(type <= 4 && type >= 1)) {
		cout 	<< "Would you like to find: " 
				<< endl << "1. a single particle energy" 
				<< endl << "2. a range" 
				<< endl << "3. reaction rate" 
				<< endl << "4. Default Preset?" << endl;
		cin >> type;
		if (!(type <= 4 && type >= 1))
		{
			cout << "error, invalid response" << endl;
		}
	}
	if (type == 4)
	{
		cout << "D-D Fusion, E=1:3keV, Ue = 27eV" << endl;
		FuseCalc.EnergyRange(2, 2, 1, 1, 1000,3000, 25);
		return 0;
	}
	while(preset != "y" && preset != "n") {
		cout << "would you like to use a particle preset? y/n" << endl;
		cin >> preset;
		if (preset != "y" && preset != "n")
		{
			cout << "error, wrong input type" << endl;
		}		
	}
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
				SvarArray.push_back(0);
				SvarArray.push_back(1);
				
			}
			else if(presetType == 2){				
				break;
			}
			else{
				cout << "error, invalid response" << endl;
			}		
		}
	}
	else if (preset == "n")
	{
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
		cout << "please enter E in eV (Energy of incoming particle, not com energy)" << endl;
		cin >> E;
		cout << "please enter Ue in eV" << endl;
		cin >> Ue;
		FuseCalc.SingleEnergy(m1, m2, Z1, Z2, E, Ue, SvarArray);
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
		FuseCalc.EnergyRange(m1, m2, Z1, Z2, Emin, Emax, Ue);
	}
	else if (type == 3)
	{
		int E;
		double Ue;
		double n1;
		double n2;
		cout << "please enter E in eV (Energy of incoming particle, not com energy)" << endl;
		cin >> E;
		cout << "please enter Ue in eV" << endl;
		cin >> Ue;
		cout << "please enter n1" << endl; //
		cin >> n1;
		cout << "please enter n2" << endl; //Paladium density 112976 n/m^3, assuming 1 deuterium per pal atom on surface
		cin >> n2;
		FuseCalc.SingleEnergy(m1, m2, Z1, Z2, E, Ue, SvarArray);

	}
	else {
		cout << "ERROR" <<endl;
		return 0;
	}	
}
#include "FusionCalculator.h"
#include "FusionCrossSectionSolver.h"
#include "Maxwellian.h"

using namespace std;

static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;
Maxwellian Maxwell;

void FusionCalculator::SingleEnergy (double m1, double m2, int Z1, int Z2, int E, double Ue, int type) {

	FCSSolver.calculateFusionCross(E, 0, m1, m2, Z1, Z2, type);
	cout << "Unshielded Sfactor: "				<< FCSSolver.S << endl;
	double fusionProb 				= FCSSolver.P;
	FCSSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, type);
	FusionCalculator::shieldedFusionProb		= FCSSolver.P;
	double mr 						= FCSSolver.mr;
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Shielded Sfactor: "	<< FCSSolver.S << endl;
	cout << "Unshielded: " 			<< fusionProb << endl;
	cout << "Shielded: " 			<< shieldedFusionProb << endl;
	cout << "mr: "					<< mr << "kg" << endl;
	cout << "Eg: "					<< FCSSolver.Eg << endl;
	cout << "E: "					<< E  << "eV" << endl;
	cout << "E/Ue: " 				<< E/Ue << endl;
	cout << "EnhancementFactor: " 	<< enhancementFactor << endl;
}

void FusionCalculator::EnergyRange (double m1, double m2, double Z1, double Z2, int Emin, int Emax, int Erange, double Ue, int type) {
	int    Esize					= Emax-Emin+1;
	double fusionProb 				[Esize];
	double shieldedFusionProb		[Esize];
	double enhancementFactors		[Esize];
	double sFactors					[Esize];
	int E = 1;
	int Eindex = 0;
	string typeString;
	switch(type){
		case 0:
			typeString = "DD";
			break;
		case 2:
			typeString = "DTr";
			break;
		case 3:
			typeString = "DHe";
			break;
	}
	string UnshieldedProbsDest 		= (string(typeString) + "RangeOutput/UnshieldedProbs.txt");
	string ShieldedProbsDest 		= (string(typeString) + "RangeOutput/ShieldedProbs.txt");
	string EnergiesDest 			= (string(typeString) + "RangeOutput/Energies.txt");
	string EnhancementFactorsDest 	= (string(typeString) + "RangeOutput/EnhancementFactors.txt");
	string SFactorsDest 			= (string(typeString) + "RangeOutput/SFactors.txt");

	
	ofstream UnshieldedProbs (UnshieldedProbsDest.c_str());
	ofstream ShieldedProbs (ShieldedProbsDest.c_str());
	ofstream Energies (EnergiesDest.c_str());
	ofstream EnhancementFactors (EnhancementFactorsDest.c_str());
	ofstream SFactors (SFactorsDest.c_str());
	int Estep = (int)((Emax-Emin)/Erange+0.5);
	for (int E = Emin; E <= Emax; E=E+Estep)
	{
		FCSSolver.calculateFusionCross(E, 0, m1, m2, Z1, Z2, type);
		fusionProb[Eindex]			=  	FCSSolver.P;
		FCSSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, type);
		shieldedFusionProb[Eindex] 	=  	FCSSolver.P;
		enhancementFactors[Eindex]	= 	shieldedFusionProb[Eindex]/fusionProb[Eindex];
		sFactors[Eindex]			= 	FCSSolver.S;
		UnshieldedProbs 			<< 	fusionProb[Eindex] 			<< endl;
		ShieldedProbs				<< 	shieldedFusionProb[Eindex] 	<< endl;
		EnhancementFactors			<< 	enhancementFactors[Eindex] 	<< endl;
		SFactors 					<<  sFactors[Eindex]			<< endl;
		Energies					<< 	FCSSolver.E 				<< endl;
		Eindex						= 	(E - Emin)/Estep;
	}
	cout << "Done, see text files in ShieldingTests/" + string(typeString) + "/RangeOutput" << endl;
}

double FusionCalculator::ReactionRate (double n1, double n2, double m1, double m2, int Z1, int Z2, int E, double Ue, int type) {	
	SingleEnergy(m1, m2, Z1, Z2, E, Ue, type);
	double vr = Maxwell.ConvertKineticEnergyToSpeed(E*eVtoJ, mr);
	double rate = n1*n2*shieldedFusionProb*vr;
	cout << "Rate = " << rate << endl; //10^14/10^16ish density for incoming beam of deuterium	
	return rate;
}
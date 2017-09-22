#include "FusionCalculator.h"
#include "FusionCrossSectionSolver.h"
#include "Maxwellian.h"

using namespace std;

static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;
Maxwellian Maxwell;

void FusionCalculator::SingleEnergy (double m1, double m2, int Z1, int Z2, int E, double Ue,
	double A1, double A2, double A3, double A4, double A5, double B1, double B2, double B3, double B4) {

	double Ej 						= E * eVtoJ;
	double Uej						= Ue *eVtoJ;

	FCSSolver.setSFactorVars(A1, A2, A3, A4, A5, B1, B2, B3, B4);
	FCSSolver.calculateFusionCross(Ej, m1, m2, Z1, Z2);
	double fusionProb 				= FCSSolver.P;
	FCSSolver.calculateFusionCross(Ej+Uej, m1, m2, Z1, Z2);
	double shieldedFusionProb		= FCSSolver.P;
	double mr 						= FCSSolver.mr;
	double Er 						= FCSSolver.Er;
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Sfactor: "				<< FCSSolver.S;
	cout << "Unshielded: " 			<< fusionProb << endl;
	cout << "Shielded: " 			<< shieldedFusionProb << endl;
	cout << "mr: "					<< mr << "kg" << endl;
	cout << "E: "					<< E  << "eV" << endl;
	cout << "Er: "					<< Er << "J" << endl;
	cout << "Er: "					<< Er/eVtoJ << "eV" << endl;
	cout << "Er/Ue: " 				<< Er/Uej << endl;
	cout << "EnhancementFactor: " 	<< enhancementFactor << endl;
}

void FusionCalculator::EnergyRange (double m1, double m2, double Z1, double Z2, int Emin, int Emax, double Ue) {
	// double Uej						= Ue *eVtoJ;
	// int    Esize					= Emax-Emin;
	// double fusionProb 				[Esize];
	// double shieldedFusionProb		[Esize];
	// double enhancementFactors		[Esize];
	// double Ej;
	// int    Eindex;
	// ofstream UnshieldedProbs ("UnshieldedProbs.txt");
	// ofstream ShieldedProbs ("ShieldedProbs.txt");
	// ofstream Energies ("Energies.txt");
	// ofstream EnhancementFactors ("EnhancementFactors.txt");
	// for (int E = Emin; E <= Emax; ++E)
	// {
	// 	Ej 							= E * eVtoJ;
	// 	Eindex						= E - Emin;
	// 	fusionProb[Eindex]			= FCSSolver.returnFusionCross(Ej, m1, m2, Z1, Z2);
	// 	shieldedFusionProb[Eindex] 	= FCSSolver.returnFusionCross(Ej+Uej, m1, m2, Z1, Z2);
	// 	enhancementFactors[Eindex]	= shieldedFusionProb[Eindex]/fusionProb[Eindex];
	// 	UnshieldedProbs 			<< fusionProb[Eindex] << endl;
	// 	ShieldedProbs				<< shieldedFusionProb[Eindex] << endl;
	// 	EnhancementFactors			<< enhancementFactors[Eindex] << endl;
	// 	Energies					<< Ej << endl;
	// }
	// cout << "Done, see text files in ShieldingTests";
}

double FusionCalculator::ReactionRate (double n1, double n2, double m1, double m2, int Z1, int Z2, double E, double Ue, double S) {	
	// SingleEnergy(m1, m2, Z1, Z2, E, Ue, A1, A2, A3, A4, A5, B1, B2, B3, B4);
	// double vr = Maxwell.ConvertKineticEnergyToSpeed(Er, mr);
	// double rate = n1*n2*shieldedFusionProb*vr;
	// cout << "Rate = " << rate << endl;
	// return rate;
}
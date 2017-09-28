#include "FusionCalculator.h"
#include "FusionCrossSectionSolver.h"
#include "Maxwellian.h"

using namespace std;

static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;
Maxwellian Maxwell;

void FusionCalculator::SingleEnergy (double m1, double m2, int Z1, int Z2, int E, double Ue, std::vector<int> SvarArray) {

	double Ej 						= E * eVtoJ;
	double Uej						= Ue *eVtoJ;
	FCSSolver.calculateFusionCross(Ej, 0, m1, m2, Z1, Z2, SvarArray);
	double fusionProb 				= FCSSolver.P;
	FCSSolver.calculateFusionCross(Ej, Uej, m1, m2, Z1, Z2, SvarArray);
	FusionCalculator::shieldedFusionProb		= FCSSolver.P;
	double mr 						= FCSSolver.mr;
	double Er 						= FCSSolver.Er;
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Sfactor: "				<< FCSSolver.S << endl;
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

double FusionCalculator::ReactionRate (double n1, double n2, double m1, double m2, int Z1, int Z2, int E, double Ue, std::vector<int> SvarArray) {	
	SingleEnergy(m1, m2, Z1, Z2, E, Ue, SvarArray);
	double vr = Maxwell.ConvertKineticEnergyToSpeed(Er, mr);
	cout << "shieldedFusionProb " << FusionCalculator::shieldedFusionProb << endl;
	double rate = n1*n2*shieldedFusionProb*vr;
	cout << "Rate = " << rate << endl; //10^14/10^16ish density for incoming beam of deuterium	
	return rate;
}
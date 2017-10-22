#include "FusionCalculator.h"
#include "FusionCrossSectionSolver.h"
#include "Maxwellian.h"

using namespace std;

static const double eVtoJ 		= 1.60218e-19;

FusionCrossSectionSolver FCSSolver;
Maxwellian Maxwell;

void FusionCalculator::SingleEnergy (double m1, double m2, int Z1, int Z2, int E, double Ue, std::vector<int> SvarArray) {

	FCSSolver.calculateFusionCross(E, 0, m1, m2, Z1, Z2, SvarArray);
	cout << "Unshielded Sfactor: "				<< FCSSolver.S << endl;
	double fusionProb 				= FCSSolver.P;
	FCSSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, SvarArray);
	FusionCalculator::shieldedFusionProb		= FCSSolver.P;
	double mr 						= FCSSolver.mr;
	double Er 						= FCSSolver.Er;
	double enhancementFactor		= shieldedFusionProb/fusionProb;
	cout << "Shielded Sfactor: "	<< FCSSolver.S << endl;
	cout << "Unshielded: " 			<< fusionProb << endl;
	cout << "Shielded: " 			<< shieldedFusionProb << endl;
	cout << "mr: "					<< mr << "kg" << endl;
	cout << "Eg: "					<< FCSSolver.Eg << endl;
	cout << "E: "					<< E  << "eV" << endl;
	cout << "Er: "					<< Er << "eV" << endl;
	cout << "Er: "					<< Er*eVtoJ << "J" << endl;
	cout << "Er/Ue: " 				<< Er/Ue << endl;
	cout << "EnhancementFactor: " 	<< enhancementFactor << endl;
}

void FusionCalculator::EnergyRange (double m1, double m2, double Z1, double Z2, int Emin, int Emax, double Ue, std::vector<int> SvarArray) {
	int    Esize					= Emax-Emin+1;
	double fusionProb 				[Esize];
	double shieldedFusionProb		[Esize];
	double enhancementFactors		[Esize];
	int E = 1;
	int Eindex;
	ofstream UnshieldedProbs ("RangeOutput/UnshieldedProbs.txt");
	ofstream ShieldedProbs ("RangeOutput/ShieldedProbs.txt");
	ofstream Energies ("RangeOutput/Energies.txt");
	ofstream EnhancementFactors ("RangeOutput/EnhancementFactors.txt");
	for (int E = Emin; E <= Emax; E++)
	{
		Eindex						= E - Emin;
		FCSSolver.calculateFusionCross(E, 0, m1, m2, Z1, Z2, SvarArray);
		fusionProb[Eindex]			=  FCSSolver.P;
		FCSSolver.calculateFusionCross(E, Ue, m1, m2, Z1, Z2, SvarArray);
		shieldedFusionProb[Eindex] 	=  FCSSolver.P;
		enhancementFactors[Eindex]	= shieldedFusionProb[Eindex]/fusionProb[Eindex];
		UnshieldedProbs 			<< fusionProb[Eindex] << endl;
		ShieldedProbs				<< shieldedFusionProb[Eindex] << endl;
		EnhancementFactors			<< enhancementFactors[Eindex] << endl;
		Energies					<< FCSSolver.Er << endl;
	}
	cout << "Done, see text files in ShieldingTests";
}

double FusionCalculator::ReactionRate (double n1, double n2, double m1, double m2, int Z1, int Z2, int E, double Ue, std::vector<int> SvarArray) {	
	SingleEnergy(m1, m2, Z1, Z2, E, Ue, SvarArray);
	double vr = Maxwell.ConvertKineticEnergyToSpeed(Er, mr);
	double rate = n1*n2*shieldedFusionProb*vr;
	cout << "Rate = " << rate << endl; //10^14/10^16ish density for incoming beam of deuterium	
	return rate;
}
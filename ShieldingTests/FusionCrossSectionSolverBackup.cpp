//Program that finds Fusion cross section for a set of particles of given relative energy 

#include "FusionCrossSectionSolver.h"

using namespace std;

const double c = 299792458;
const double pi = 3.14159265359;
const double eVtoJ = 1.60218e-19;
const double alpha = 1/137.035999;

void FusionCrossSectionSolver::calculateReducedMass(double m1, double m2) {
	FusionCrossSectionSolver::mr = (m1 * m2)/(m1 + m2);
	// cout << "mr = " << mr << " m1 = " << m1 << " m2 = " << m2 << endl;
}

void FusionCrossSectionSolver::calculateReducedEnergy(double m1, double E){
	FusionCrossSectionSolver::Er = mr/m1*E;
}
void FusionCrossSectionSolver::calculateGamowEnergy(double mr, int Za, int Zb) {
	double subVar1 = pi * alpha * Za * Zb; // /137 due to alpha
	FusionCrossSectionSolver::Eg = 2 * mr * pow(c, 2) * pow(subVar1,2)/eVtoJ;
}
void FusionCrossSectionSolver::setSFactorVars(int fusionType) {
	switch(fusionType){
		case 0:
			//D(d,p)Tr
			FusionCrossSectionSolver::A1 = 5.5576e4;
			FusionCrossSectionSolver::A2 = 2.1054e2;
			FusionCrossSectionSolver::A3 = -3.2638e-2;
			FusionCrossSectionSolver::A4 = 1.4987e-6;
			FusionCrossSectionSolver::A5 = 1.8181e-10;
			FusionCrossSectionSolver::B1 = 0;
			FusionCrossSectionSolver::B2 = 0;
			FusionCrossSectionSolver::B3 = 0;
			FusionCrossSectionSolver::B4 = 0;
			break;
		case 1:
			//D(d,n)3He
			FusionCrossSectionSolver::A1 = 5.3701e4;
			FusionCrossSectionSolver::A2 = 3.3027e2;
			FusionCrossSectionSolver::A3 = -1.2706e-1;
			FusionCrossSectionSolver::A4 = 2.9327e-5;
			FusionCrossSectionSolver::A5 = -2.5151e-9;
			FusionCrossSectionSolver::B1 = 0;
			FusionCrossSectionSolver::B2 = 0;
			FusionCrossSectionSolver::B3 = 0;
			FusionCrossSectionSolver::B4 = 0;
			break;
		case 2:
			// D(t,n)alpha
			FusionCrossSectionSolver::A1 = 6.927e4;
			FusionCrossSectionSolver::A2 = 7.454e8;
			FusionCrossSectionSolver::A3 = 2.05E6;
			FusionCrossSectionSolver::A4 = 5.2002e4;
			FusionCrossSectionSolver::A5 = 0;
			FusionCrossSectionSolver::B1 = 6.38e1;
			FusionCrossSectionSolver::B2 = -9.95e-1;
			FusionCrossSectionSolver::B3 = 6.981e-5;
			FusionCrossSectionSolver::B4 = 1.728e-4;
			break;
		case 3:
			//D(3He,n)alpha
			FusionCrossSectionSolver::A1 = 5.7501e16;
			FusionCrossSectionSolver::A2 = 2.5226e3;
			FusionCrossSectionSolver::A3 = 4.5566e1;
			FusionCrossSectionSolver::A4 = 0;
			FusionCrossSectionSolver::A5 = 0;
			FusionCrossSectionSolver::B1 = -3.1995e-3;
			FusionCrossSectionSolver::B2 = -8.5530e-6;
			FusionCrossSectionSolver::B3 = 5.9014e-8;
			FusionCrossSectionSolver::B4 = 0;
		break;
	}
}
void FusionCrossSectionSolver::calculateSFactor(double Ue, std::vector<int> SvarArray) {
	FusionCrossSectionSolver::S = 0;
	for (int i = 0; i < (SvarArray.size()); ++i)
	{
		setSFactorVars(SvarArray[i]);
		double ErkeV = (Er + Ue)/(1000);
		double Snumerator = (A1+ErkeV*(A2+ErkeV*(A3+ErkeV*(A4 + ErkeV *A5))));
		double Sdenominator = (1+ErkeV*(B1+ErkeV*(B2+ErkeV*(B3 + ErkeV * B4))));
		FusionCrossSectionSolver::S += (Snumerator/Sdenominator);
		// cout << "A1= " << A1 << endl;
		// cout << "Snumerator = " << Snumerator << endl;
		// cout << "SinEv = " << S << endl;
	}	
	FusionCrossSectionSolver::S = (S/SvarArray.size())*1000;
	// cout << "finalS = " << S << endl;
}
	
void FusionCrossSectionSolver::calculateFusionCross(double E, double Ue, double m1, double m2, int Za, int Zb, std::vector<int> SvarArray) {
	calculateReducedMass(m1, m2);
	calculateGamowEnergy(mr, Za, Zb);
	calculateReducedEnergy(m1,E);
	calculateSFactor(Ue, SvarArray);
	double SofEonE =S/(Er+Ue);
	// double SofEonE = 1;
	FusionCrossSectionSolver::P = SofEonE*exp(-sqrt(Eg/(Er+Ue))); //*10^-28
	//cout << "T = " << T << endl;
}
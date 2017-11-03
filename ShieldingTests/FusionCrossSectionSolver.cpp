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
			FusionCrossSectionSolver::A1 = 5.7501e6;
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
void FusionCrossSectionSolver::calculateSFactor(int type) {
	FusionCrossSectionSolver::S = 0;
	bool SfactorDone = false;
	while(!SfactorDone){
		setSFactorVars(type);
		double EkeV = (E)/(1000);
		double Snumerator = (A1+EkeV*(A2+EkeV*(A3+EkeV*(A4 + EkeV *A5))));
		double Sdenominator = (1+EkeV*(B1+EkeV*(B2+EkeV*(B3 + EkeV * B4))));
		FusionCrossSectionSolver::S += (Snumerator/Sdenominator)*1000;
		// cout << "A1= " << A1 << endl;
		// cout << "Snumerator = " << Snumerator << endl;
		// cout << "SinEv = " << S << endl;
		if (type == 0)
		{
			type = 1;
		}
		else{
			SfactorDone = true;
		}
	}	
	if (type == 1)
	{
		FusionCrossSectionSolver::S = (FusionCrossSectionSolver::S/2);
	}
	// cout << "finalS = " << S << endl;
}

	
void FusionCrossSectionSolver::calculateFusionCross(double E, double Ue, double m1, double m2, int Za, int Zb, int type) {
	FusionCrossSectionSolver::E = E;
	calculateReducedMass(m1, m2);
	calculateGamowEnergy(mr, Za, Zb);
	calculateSFactor(type);
	double SofEonE =S/sqrt(E*(E+Ue));
	// double SofEonE = 1;
	FusionCrossSectionSolver::P = SofEonE*exp(-sqrt(Eg/(E+Ue))) *1e-31; //*10^-28
	//cout << "T = " << T << endl;
}
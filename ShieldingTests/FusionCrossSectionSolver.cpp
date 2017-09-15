//Program that finds Fusion cross section for a set of particles of given relative energy 

#include "FusionCrossSectionSolver.h"

using namespace std;

const double c = 299792458;
const double pi = 3.14159265359;
const double eVtoJ = 1.60218e-19;
const double alpha = 1/137.035999;

double FusionCrossSectionSolver::returnReducedMass(double m1, double m2) {
	double mr = (m1 * m2)/(m1 + m2);
	// cout << "mr = " << mr << " m1 = " << m1 << " m2 = " << m2 << endl;
	return mr;
}

double FusionCrossSectionSolver::returnGamowEnergy(double mr, int Za, int Zb) {
	double subVar1 = pi * alpha * Za * Zb; // /137 due to alpha
	double Eg = 2 * mr * pow(c, 2) * pow(subVar1,2);
	// cout << "Eg = " << Eg << endl;
	return Eg;
}
	
double FusionCrossSectionSolver::returnFusionCross(double E, double m1, double m2, int Za, int Zb) {
	double mr 	= returnReducedMass(m1, m2);
	double Eg 	= returnGamowEnergy(mr, Za, Zb);
	double P 	= exp(-sqrt(Eg/E));
	//cout << "T = " << T << endl;
	return P;
}

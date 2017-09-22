//Program that finds Fusion cross section for a set of particles of given relative energy 

#include "FusionCrossSectionSolver.h"

using namespace std;

const double c = 299792458;
const double pi = 3.14159265359;
const double eVtoJ = 1.60218e-19;
const double alpha = 1/137.035999;

void FusionCrossSectionSolver::calculateReducedMass(double m1, double m2) {
	double mr = (m1 * m2)/(m1 + m2);
	// cout << "mr = " << mr << " m1 = " << m1 << " m2 = " << m2 << endl;
}

void FusionCrossSectionSolver::calculateReducedEnergy(double m1, double mr, double E){
	double Er = mr/m1*E;
}
void FusionCrossSectionSolver::FusionCrossSectionSolver::calculateGamowEnergy(double mr, int Za, int Zb) {
	double subVar1 = pi * alpha * Za * Zb; // /137 due to alpha
	double Eg = 2 * mr * pow(c, 2) * pow(subVar1,2);
	// cout << "Eg = " << Eg << endl;
}
void FusionCrossSectionSolver::setSFactorVars(double A1par, double A2par, double A3par, double A4par, double A5par, double B1par, double B2par, double B3par, double B4par){
	double A1 = A1par;
	double A2 = A2par;
	double A3 = A3par;
	double A4 = A4par;
	double A5 = A5par;
	double B1 = B1par;
	double B2 = B2par;
	double B3 = B3par;
	double B4 = B4par;
}
void FusionCrossSectionSolver::calculateSFactor(double Eg, double E) {
	double S = (A1+E*(A2+E*(A3+E*(A4 + E *A5)))/(1+E*(B1+E*(B2+E*(B3 + E * B4)))));
}
	
void FusionCrossSectionSolver::calculateFusionCross(double E, double m1, double m2, int Za, int Zb) {
	calculateReducedMass(m1, m2);
	calculateGamowEnergy(mr, Za, Zb);
	calculateReducedEnergy(m1, mr, E);
	calculateSFactor(Eg,Er);
	SofEonE = S/E;
	double P 	= SofEonE*exp(-sqrt(Eg/Er));
	//cout << "T = " << T << endl;
}
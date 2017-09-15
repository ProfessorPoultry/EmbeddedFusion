#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "Maxwellian.h"
#include "FusionCrossSectionSolver.h"

using namespace std;

static const double k 		 	= 1.38064852e-23;
static const double m_proton 	= 1.6726219e-27;
static const double eVtoJ 		= 1.60218e-19;
void main () {
	cout << "m_proton = " << m_proton << endl;
	double m1 = 2 * m_proton;
	double m2 = m1;
	int Za = 1;
	int	Zb = 1;
	double E = 1000* eVtoJ;
	FusionCrossSectionSolver FCSSolver;
	double FusionProb = FCSSolver.returnFusionCross(E, m1, m2, Za, Zb);
	cout << FusionProb;
}
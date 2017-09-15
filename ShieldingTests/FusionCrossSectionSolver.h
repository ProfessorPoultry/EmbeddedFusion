#ifndef FUSIONCROSSSECTIONSOLVER_H
#define FUSIONCROSSSECTIONSOLVER_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

class FusionCrossSectionSolver {
	public:
		double T;
		double mr;
		double Eg;
		double returnReducedMass(double, double);
		double returnReducedEnergy(double, double, double);
		double returnGamowEnergy(double, int, int);
		double returnFusionCross(double, double, double, int, int);
};

#endif
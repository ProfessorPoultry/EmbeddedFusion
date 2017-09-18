#ifndef FUSIONCALCULATOR_H
#define FUSIONCALCULATOR_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

class FusionCalculator {
	double Ej;
	double Uej;
	double fusionProb;
	double shieldedFusionProb;
	double mr;
	double Er;
	double enhancementFactor;
	double vr;
	double rate;

	public:
		void SingleEnergy (double , double , double , double , double,  double);
		void EnergyRange (double , double , double , double , int, int, double);
 		double ReactionRate (double , double , double , double , double , double , double, double);
};

#endif
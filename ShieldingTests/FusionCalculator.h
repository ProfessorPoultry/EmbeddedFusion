#ifndef FUSIONCALCULATOR_H
#define FUSIONCALCULATOR_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>

class FusionCalculator {
	

	public:
		double Ej;
		double Uej;
		double fusionProb;
		double shieldedFusionProb;
		double mr;
		double Er;
		double enhancementFactor;
		double vr;
		double rate;
		void SingleEnergy (double , double , int , int , int,  double, int);
		void EnergyRange (double , double , double , double , int, int, int, double, int);
 		double ReactionRate (double , double , double , double , int , int , int, double, int);
};

#endif
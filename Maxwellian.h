#ifndef MAXWELLIAN_H
#define MAXWELLIAN_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

class Maxwellian {
	double T;
	double m;
	double speeds[10000];
	double probabilities[10000];
	public:
		double GetProbability(int);
		double ConvertSpeedToKineticEnergy(double, double);
		double ConvertKineticEnergyToSpeed(double, double);
		void   PopulateMaxwellian(double, double, double, double, double);
		double GetPlasmaTemp();
		double GetSpeed(int);
		double CalculateProbability(int, double, double);
};

#endif
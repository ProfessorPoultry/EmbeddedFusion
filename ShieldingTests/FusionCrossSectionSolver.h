#ifndef FUSIONCROSSSECTIONSOLVER_H
#define FUSIONCROSSSECTIONSOLVER_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>

class FusionCrossSectionSolver {
	public:

		//S Factor variables
		double A1;
		double A2;
		double A3;
		double A4;
		double A5;
		double B1;
		double B2;
		double B3;
		double B4;

		//dependent variables
		double T;
		double mr;
		double Er;
		double Eg;
		double S;
		double SofEonE;
		double P;

		void calculateReducedMass(double, double);
		void calculateReducedEnergy(double, double, double);
		void calculateGamowEnergy(double, int, int);
		void setSFactorVars(double, double, double, double, double, double, double, double, double);
		void calculateSFactor(double, double);
		void calculateFusionCross(double, double, double, int, int);
};

#endif
#ifndef FUSIONCROSSSECTIONSOLVER_H
#define FUSIONCROSSSECTIONSOLVER_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>

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
		void calculateReducedEnergy(double, double);
		void calculateGamowEnergy(double, int, int);
		void setSFactorVars(int);
		void calculateSFactor(double, std::vector<int>);
		void calculateFusionCross(double,double, double, double, int, int, std::vector<int>);
};

#endif
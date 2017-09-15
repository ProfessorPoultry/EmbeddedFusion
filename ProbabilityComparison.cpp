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
	FusionCrossSectionSolver FCSSolver;
	Maxwellian Plasma1;
	double T				= 100000; //Kelvin
	double m1				= 2 * m_proton;
	double m2 				= m1;
	double minEnergy		= 150 * eVtoJ;
	double maxEnergy		= 1500 * eVtoJ;
	double energyInterval	= 0.01 * eVtoJ;
	double speedInterval	= Plasma1.ConvertKineticEnergyToSpeed(energyInterval,m1);
	double min_speed		= Plasma1.ConvertKineticEnergyToSpeed(minEnergy,m1);
	double max_speed		= Plasma1.ConvertKineticEnergyToSpeed(maxEnergy,m1); //m/s
	int Za = 1;
	int	Zb = 1;

	// Plasma1.PopulateMaxwellian(T, m, min_speed, speedInterval, max_speed);
	ofstream EnergyOutput ("EnergyOutput.txt");
	ofstream PlasmaProbOutput ("PlasmaProbOutput.txt");
	ofstream FusionProbOutput ("FusionProbOutput.txt");
	ofstream CombinedProbOutput ("CombinedProbOutput.txt");
	if (EnergyOutput.is_open() && PlasmaProbOutput.is_open() && FusionProbOutput.is_open() && CombinedProbOutput.is_open()) {
  		double v = min_speed;
  		double E;
  		double FusionProb;
  		double plasmaProb;
  		double combinedProbability;

  		int index = 0;
		while(v <= max_speed) {
			E 					= Plasma1.ConvertSpeedToKineticEnergy(v, m1);
			FusionProb 		 	= FCSSolver.returnFusionCross(E, m1, m2, Za, Zb);
			plasmaProb 			= Plasma1.CalculateProbability(v,T,m1);
			combinedProbability = FusionProb * plasmaProb; //Because cross section is just tunneling probability
			
			EnergyOutput	 	<< E/eVtoJ	   			<< endl;
			PlasmaProbOutput 	<< plasmaProb 			<< endl;
			FusionProbOutput 	<< FusionProb 			<< endl;
			CombinedProbOutput 	<< combinedProbability 	<< endl;
			
			v+=speedInterval;
		}
	}	
	EnergyOutput.close();
	PlasmaProbOutput.close();
	FusionProbOutput.close();
	CombinedProbOutput.close();
	return 0;
}
	/*
		Gamow factor
		P = probability of overcoming/tunneling through coulomb barrier

		P(E) = e ^ - ((Eg/E)^0.5)
		Eg = 2 * mr * c^2 *(pi*alpha*Za*Zb)^2

		mr = reduced mass = m1*m2/(m1+m2)

		Za and Zb = atomic numbers

		c = speed of light

		alpha = fine structure constant
	*/
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "Maxwellian.h"

using namespace std;

static const double k 		 	= 1.38064852e-23;
static const double m_proton 	= 1.6726219*10e-27;
static const double eVtoJ 		= 1.60218e-19;

double ConvertKineticEnergyToSpeed(double E, double m) {
	double v = sqrt(2*E/m);
	return v;
}

double ConvertSpeedToKineticEnergy(double v, double m) {
	double E = pow(v,2)*m/2;
	return E;
}

int main() {

	Maxwellian Plasma1;

	double T 	= 100000; //Kelvin
	double m 	= 2* m_proton;
	double minEnergy = 0*eVtoJ;
	double maxEnergy = 6000*eVtoJ;
	double energyInterval = 0.1*eVtoJ;
	double speedInterval  = ConvertKineticEnergyToSpeed(energyInterval,m);
	double min_speed 	= ConvertKineticEnergyToSpeed(minEnergy,m);
	double max_speed 	= ConvertKineticEnergyToSpeed(maxEnergy,m); //m/s

	//Plasma1.PopulateMaxwellian(T, m, min_speed, speedInterval, max_speed);
	ofstream output ("output.txt");
  	if (output.is_open()) {
  		double v = min_speed;
  		cout << v << endl;
		while(v <= max_speed) {
			output << Plasma1.GetProbability(v,T,m) << endl;
			v+=speedInterval;
		}
		cout << speedInterval << endl << v << endl;
	}	
	output.close();
	double onekEVSpeed = ConvertKineticEnergyToSpeed(1000*eVtoJ, m);
	double onekEVProb = Plasma1.GetProbability(onekEVProb,T,m);
	cout << onekEVSpeed << endl;
	cout << onekEVProb;
	return 0;
}


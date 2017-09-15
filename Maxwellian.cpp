//Maxwellian probability
#include "Maxwellian.h"

using namespace std;

const double k = 1.38064852e-23;
const double pi = 3.14159265359;
const double eVtoJ = 1.60218e-19;

double Maxwellian::ConvertKineticEnergyToSpeed(double E, double m) {
	double v = sqrt(2*E/m);
	return v;
}

double Maxwellian::ConvertSpeedToKineticEnergy(double v, double m) {
	double E = pow(v,2)*m/2;
	return E;
}

double Maxwellian::GetProbability(int i){
	// double probabilities;
	// return probabilities[i];
}
void Maxwellian::PopulateMaxwellian(double T, double m, double min_speed, double speedInterval, double max_speed) {
	double current_speed = 0;
	double v = min_speed;
	double x = 0;
	int index = 0;
	while(v <= max_speed) {
		// probabilities[v] = pow(m/(2*pi*k*T), (3/2)) * 4 * pi * pow(v,2) * exp(-m*pow(v,2)/(2*k*T));
		x = v/sqrt(2*k*T/m);
		probabilities[index] = CalculateProbability(x, T, m);
		v+=speedInterval;
		index++;
	}
}

double Maxwellian::GetPlasmaTemp() {
	return T;
}

double Maxwellian::GetSpeed(int i) {
	return speeds[i];
}

double Maxwellian::CalculateProbability(int v, double T, double m) {
	//return 4/sqrt(pi)*(pow(x,2))*exp(-pow(x,2)); //dimensionless
	return sqrt(pow(m/(2*pi*k*T),3))*4*pi*pow(v,2)*exp(-m*pow(v,2)/(2*k*T));
}

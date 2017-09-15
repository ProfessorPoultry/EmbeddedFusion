// This is a Poisson solver in cylindrical symmetry for a two ring cathode 

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main() {

// Potiential is the 2D cylindrical array of potentials
int i, j, k, mesh_x=100, mesh_y=50, radius, ring1_z, ring2_z;
double potential[101][51], boundary[101][51], old_potential[101][51], density[101][51], difference, old_difference, tolerence, max_difference, temp;                 //Temporary density declaration
   
//denisty = 1.6e-19 * (ni + ne) / 8.85e-12                               //Where ni == ne

radius = 1;
ring1_z = 40;
ring2_z = 60;
tolerence = 0.0001; // Tolerence in volts
max_difference = tolerence+1; // initial differences in maximum change between old and new potentials

// Set the boundary conditions
// This is two ring grid

// zero entire array
for (i = 0; i<=mesh_x; i++){
    for (j = 0; j<=mesh_y; j++){    
        potential[i][j]=0;
        density[i][j]=0;
        old_potential[i][j]=0;
        }
}


// The potential V=-1000V on the two rings. Set these as boundaries.
potential[ring1_z][radius] = -1000;
boundary[ring1_z][radius] = 1;
potential[ring2_z][radius] = -1000;
boundary[ring2_z][radius] = 1;

while (max_difference > tolerence){

//The start of the algorithm to solve Poisson's equation. The outer boundaries are left alone so the loops do not start at 0 and end on 50 or 100.
      for (i=1; i<mesh_x; i++){	
          for (j=1; j<mesh_y; j++){	
		
// Check if you are at the rings. If so then skip the algorithm		
                if (boundary[i][j] !=1 ){
		
// The algorithm of averaging nearest neighbours in cylindrical coordinates
                    potential[i][j] = potential[i+1][j]*(1+1/(2*i)) + potential[i-1][j]*(1-1/(2*i)) + potential[i][j+1] + potential[i][j-1] ;
                    potential[i][j]= potential[i][j]/4;
                    }	
                 }
          }

// Now assign potentials along the axis using the Neumman condition.
   for (i=1; i<mesh_x; i++){
       potential[i][0] = potential[i][1];	
   }

// Obtain the total absolute difference between the old an new potentials. Then transfer the current potential to the old potential.
// Obtain the maximum difference between old and new to compare it to the tolerence at the start of the main loop
max_difference = 0;
for (i=0; i<=mesh_x; i++){
	for (j=0; j<=mesh_y; j++){
		temp = potential[i][j] - old_potential[i][j];
		difference = sqrt(temp*temp);
		if (difference > max_difference){
			max_difference = difference;
			
		}
	}
}

// Reset the difference variable for the next time around in the loop.
// Transfer all the current potentials into the old_potential so that they can be compared to 
// the updated potentials next time around in the loop

difference = 0;
for (i=0; i<=mesh_x; i++){
	for (j=0; j<=mesh_y; j++){
		old_potential[i][j] = potential[i][j];
	}
}
}


// Output potential along the axis to a file
 ofstream myfile;
  myfile.open ("poisson.csv");
  for (i=0; i<=mesh_x; i++){  
  myfile << potential[i][0] << endl;
 
  	 }

   myfile.close();
   
          system("pause");

}

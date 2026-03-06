#include <stdio.h>
#include <math.h>


// number of points in the potential
const int N=100;

// the things needed to define the potential
double r[N]; // distances between particles
double U[N]; // potential at distances defiend in r
double F[N]; // force at distances defined in r

double r_cut = 5;

double sigma = 1;
double ebsilon = 1;



void lenard_jones(){
    double base_c = sigma/r_cut;
    double e_cut = ebsilon*4*  ( pow( base_c, 12) - pow(base_c,6) );


    // 0 has an undefined potential
    for(int n=1; n<N; n++){
        // generate the 3 where we reside
        r[n] = (double)n / N * r_cut;

        // generating the base
        double base = sigma / r[n];

        // defining the potential
        U[n] = ebsilon*4*  ( pow( base, 12) - pow(base,6) ) - e_cut;

        if(n>1){
            F[n] = - (U[n-1]-U[n]) / (r[n-1]-r[n]);
        }
        else{
            F[n]=pow(10,10); //potential decreases very fast at 0 distance meaning we want a very high positive number
        }
    }

}



void main(){



    
}
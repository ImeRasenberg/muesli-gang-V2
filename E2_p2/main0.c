#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 256

const char* init_filename = "fcc.xyz";

/* Simulation variables */
double dt = 1E-3; 
#define M 5000
double time_array[M];

/* System variables*/
int n_particles; // first three are imported from fcc.xyz
double radius;
double box[NDIM];

double beta = 1.0;
double mass = 1.0;
double sigma;
double density = 0.5;

double E[M][2]; // kinetic and potential energy 
double r[N][NDIM];
double v[N][NDIM];
double F[N][NDIM];

/*potential variables*/
double epsilon = 1.0;
double r_cut;
double e_cut = 0.0;
double energy;




void read_data(void){
   FILE* coords;
   coords = fopen(init_filename, "r");

   fscanf(coords, "%i\n", &n_particles); // read in number of particles

   for (int i = 0; i < NDIM; i++) {// read in dimension of box 
    double L, R;
    fscanf(coords, "%lf %lf", &L, &R); 
    box[i] = R-L;
    }

    for (int j = 0; j < n_particles; j++) { // read in partilces pos
        if (j == 0){ // use the radisu of the first line (assumed to be the same for all)
         fscanf(coords, "%lf %lf %lf %lf ", &r[j][0], &r[j][1], &r[j][2], &radius);   
        }
        else fscanf(coords, "%lf %lf %lf %*lf ", &r[j][0], &r[j][1], &r[j][2]);
    }
    fclose(coords);
    /* print out some of the read data as a sanity check*/
    printf("number of particls: %d\t ,box size: %lf\t ,radius: %lf\n " , n_particles, box[0],radius ); 
}

void set_density(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

void gaussian_numbers(double *z0, double *z1) {
    double U = dsfmt_genrand();
    double V = dsfmt_genrand();

    double r = sqrt(-2*log(U));
    double theta = 2* M_PI * V;

    *z0 = r * cos(theta);
    *z1 = r * sin(theta);

}
void initialize_velocity(void){
    
    /*intilize the box muller gausian velociut*/
    double X, Y;
    int pairs  = N*NDIM/2; // N is equal in our case
    int idx = 0;
    for (int i = 0; i<pairs; i++){
        gaussian_numbers(&X, &Y);
        v[idx/NDIM][idx %NDIM]= X* sigma;
        idx++;
        v[idx/NDIM][idx %NDIM]= Y *sigma;
        idx++;
    }

    /* setting total momentum to zero*/
    double v_avg[NDIM] = {0.0};
    for (int n = 0; n<N; n++){
        for (int d = 0; d<NDIM; d++){
            v_avg[d] += v[n][d];
        }
    }
    for(int d = 0; d<NDIM; d++){
        v_avg[d] /= N;
    } 

    for (int n = 0; n<N; n++){
        for (int d = 0; d<NDIM; d++){
            v[n][d] -= v_avg[d];
        }
    }
}

double calculate_potential(double r2) {

    if (r2 >= r_cut*r_cut) return 0.0; 

    double inv_r2 = 1.0 / r2;
    double sr2 = sigma*sigma * inv_r2;
    double sr6 = sr2*sr2*sr2;
    double sr12 = sr6 * sr6;
    
    return (4*epsilon) * ( sr12 - sr6) + epsilon;
}

double calculate_force_over_r(double r2) {
    if (r2 >= r_cut * r_cut) return 0.0; 

    double inv_r2 = 1.0 / r2;                 
    double sr2 = (sigma * sigma) * inv_r2;    
    double sr6 = sr2 * sr2 * sr2;             
    double sr12 = sr6 * sr6;                  
    
    return 24.0 * epsilon * sr2 * (2.0 * sr12 - sr6);
}

void calc_forces(){
    for(int i=0; i<n_particles; i++) 
        for(int k=0; k<NDIM; k++) F[i][k] = 0.0;

    energy = 0;
    double r_cut2 = r_cut * r_cut;

    for(int i = 0; i < n_particles; i++){
        for(int j = 0; j < i; j++){
            double r_ij[NDIM];
            double r2 = 0;
            

            for(int k=0; k < NDIM; k++) {
                r_ij[k] = r[i][k] - r[j][k];
  
                if (r_ij[k] >  0.5 * box[k]) r_ij[k] -= box[k];
                if (r_ij[k] < -0.5 * box[k]) r_ij[k] += box[k];
                r2 += r_ij[k] * r_ij[k];
            }


            if(r2 < r_cut2){
                double f_over_r = calculate_force_over_r(r2); 
                
                for(int k=0; k < NDIM; k++) {
                    double force_k = f_over_r * r_ij[k];
                    F[i][k] += force_k;
                    F[j][k] -= force_k;
                    // printf("force calculated i: %e\n",force_k);
                }
                energy += calculate_potential(r2);
            }
            
        }
    }
}

int main(void){
    sigma = sqrt(pow(beta*mass, -1));
    r_cut = pow(2.0, 1.0/6.0);
    read_data();
    if(n_particles == 0){
            printf("Error: Number of particles, n_particles = 0.\n");
            return 1;
    }
    set_density();
    initialize_velocity();
    calc_forces(); 
    
    for(int t = 0; t<M ; t++){
        time_array[t]= t*dt;
    }

    return 0;
}

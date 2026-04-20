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
#define dt 1E-4
#define t_max 2
#define M (int)(t_max/dt)
double time_array[M];

/* System variables*/
int n_particles; // first three are imported from fcc.xyz
double sigma;
double box[NDIM];

double beta = 1.0/3.0;
double mass = 1.0;

double density = 0.5;
double gamma = 1;
double D;

double E[M]; // kinetic and potential energy 
double r[N][NDIM];
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
         fscanf(coords, "%lf %lf %lf %lf ", &r[j][0], &r[j][1], &r[j][2], &sigma);   
        }
        else fscanf(coords, "%lf %lf %lf %*lf ", &r[j][0], &r[j][1], &r[j][2]);
    }
    fclose(coords);
    /* print out some of the read data as a sanity check*/
    printf("number of particls: %d\t ,box size: %lf\t ,radius: %lf\n " , n_particles, box[0],sigma ); 
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

double  gaussian_number(void) {
    double U = dsfmt_genrand();
    double V = dsfmt_genrand();

    double r = sqrt(-2*log(U));
    double theta = 2* M_PI * V;

    return r * cos(theta);

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

void write_data(void){
    char filename[100];
        sprintf(filename, "data/energy_vs_time_Temp_%.5f.txt", 1/beta);

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# Time    Energy\n");
        for (int t = 0; t < M; t++) {
            fprintf(fp, "%lf\t%lf\t%lf\n", time_array[t], E[t],1/beta);
        }
        fclose(fp);
}

void write_coords(double step){
    char buffer[128];
    sprintf(buffer, "data/coords%03lf.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", sigma);
    }
    fclose(fp);
}

void brownian (void){
    for (int n = 0; n<n_particles; n++){
        for (int d = 0; d <NDIM; d++){
            double R = gaussian_number();
            r[n][d] += sqrt(2*D) * R * sqrt(dt) + F[n][d]*dt/gamma;

            if (r[n][d] < 0) r[n][d] += box[d];
            if (r[n][d] >= box[d]) r[n][d] -= box[d];

        }
    }
}

int main(void){
    dsfmt_seed(time(NULL));
    read_data();
    if(n_particles == 0){
            printf("Error: Number of particles, n_particles = 0.\n");
            return 1;
    }
    
    D =  1/(beta * gamma);
    r_cut = pow(2.0, 1.0/6.0) * sigma;
    
    set_density();
    
    for(int t = 0; t<M ; t++){
        calc_forces();
        brownian(); 
        time_array[t]= t*dt;
        E[t]= energy;  
    }

    printf("simulation done");
    write_data();
    write_coords(1/beta);
    return 0;
}

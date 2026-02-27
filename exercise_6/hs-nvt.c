#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 10000;
const int output_steps = 1000;
const double packing_fraction = 0.35;
const double diameter = 1.0;
const double delta = 0.05;
const char* init_filename = "fcc.xyz";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];



/* Functions */
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
    fscanf(coords, "%lf %lf %lf %*lf ", &r[j][0], &r[j][1], &r[j][2]);
    }
    fclose(coords);
    /* print out some of the read data as a sanity check*/
    printf("number of particls %d\n", n_particles);
    printf("box size %lf\n", box[0]);
    for (int i = 0; i<3; i++){
        printf("partilce %d\t ,positions: %lf\t%lf\t%lf\n",i, r[i][0], r[i][1], r[i][2] );
    }

}

int move_particle(void){
/*--------- Your code goes here -----------*/
    int rand_particle = (int)(n_particles * dsfmt_genrand()); // choose random particle
    double trail_pos [NDIM]; // trail vector 
   
    /*move random particle */
    for (int i= 0; i<NDIM; i++){
        double rand_delta = 2 * delta* dsfmt_genrand() - delta; //generate rand [-delta, delta]
        trail_pos[i] = r[rand_particle][i]+ rand_delta; // trail move particle

        if (trail_pos[i]< 0) trail_pos[i] += box[i]; // enforce periodicity
        else if(trail_pos[i]> box[i]) trail_pos[i] -= box[i];
    }
    
    /*checking overlaps*/
    for (int n = 0; n< n_particles; n++){
    
        if (n!= rand_particle){ // skip the particle we are moving
            double distance_2 = 0;
            for (int d = 0; d< NDIM; d++){ 
                double dx = r[n][d] - trail_pos[d];
                dx = dx - box[d] * round(dx / box[d]); // nearest image
                distance_2 += dx * dx;
            }
             if (distance_2 < diameter*diameter) return 0;
        }
    }
    /*updating r vector if the trail move is accepted */
    for (int d= 0; d<NDIM ;d++)r[rand_particle][d] = trail_pos[d];
    return 1;
}

void write_data(void){
    char buffer[128];
    sprintf(buffer, "final.dat");
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

int main(int argc, char* argv[]){


    assert(packing_fraction > 0.0 && packing_fraction < 1.0);
    assert(diameter > 0.0);
    assert(delta > 0.0);

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_packing_fraction();

    dsfmt_seed(time(NULL));

    int accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }

        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %lf.\n", step, (double)accepted / (n_particles * output_steps));
            accepted = 0;
            
        }

        if ((step+1)% mc_steps== 0)write_data(); 
    }

    return 0;
}

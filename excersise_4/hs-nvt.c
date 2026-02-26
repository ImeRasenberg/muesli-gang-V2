#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "../downloads/mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3

/* Initialization variables */
const int mc_steps = 1000;
const int output_steps = 100;
const double packing_fraction = 0.55;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "fcc.xyz";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double (*r)[3];
double *size;
double box[NDIM];

// randome movement
double dr[3];
int n;

double dummy;


void read_data(void){

    FILE *read_cords;
    read_cords = fopen(init_filename, "r");

    fscanf(read_cords, "%d\n", &n_particles); 

    r = malloc(n_particles * sizeof * r); 
    size = malloc(n_particles * sizeof * size); 

    for(int i = 0; i<NDIM; i++){
        fscanf(read_cords, "%lf %lf", &dummy, &box[i]);

    }

    for(int i = 0; i<n_particles; i++){
        fscanf(read_cords, "%lf %lf %lf %lf", &r[i][0], &r[i][1], &r[i][2], &size[i]);
    }

    fclose(read_cords);
}


int check_particle_overlap(int n){
    double *p = r[n];

    for(int i=0;i<n_particles;i++){
        if(i==n){
            continue;
        }

        double *p_c = r[i];
        double distance_squared = 0;

        for(int j=0; j<3; j++){
            double difference = p[j] + dr[j] - (p_c[j]);
            if(difference>0.5*box[j]){
                difference -= box[j];
            }
            else if(difference<-0.5*box[j]){
                difference += box[j];
            }
            distance_squared += difference*difference;
        }

        double sum_raduss = 0.5 * (size[i] + size[i]);

        

        if (distance_squared < sum_raduss*sum_raduss){
            // printf("%lf\t < \t %lf\n",distance_squared,sum_raduss*sum_raduss);
            return 1;
        }
    }
    // printf("accept\n");
    return 0;
}


int move_particle(void){
    n = floor(dsfmt_genrand()*n_particles);
double dummy;
    for(int i=0;i<3;i++){
        dr[i] = (dsfmt_genrand()-0.5) + 0.00001;
    }
    double delta_l=(dsfmt_genrand()-0.5)*2*delta + 0.00001;

    double length = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    for(int i=0;i<3;i++){
        dr[i] *= delta_l/length;
    }

    int disp = check_particle_overlap(n);
    
    if(disp ==1){
        return 0;
    }
    else if (disp==0){
        for(int i=0;i<3;i++){
            r[n][i] += dr[i];

            if(r[n][i]<0){
                r[n][i] +=box[i];
            }
            if(r[n][i]>box[i]){
                r[n][i] -=box[i];
            }
        }
        return 1;
    }

}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "./data/coords_step%07d.dat", step);
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
        // printf("%lf\n",scale_factor);
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

int main(int argc, char* argv[]){
    read_data();
    
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
            write_data(step);
        }
    }
    printf("done");
    return 0;
}

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "../downloads/mt19937.h"

#include <sys/stat.h>
#include <sys/types.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512


/* Initialization variables */
const int    mc_steps      = 5000;
const int    output_steps  = 100; 
// const double density       = 0.5; 
const double delta         = 0.1; 
const double r_cut         = 2.5; 
// const double beta          = 1;
const char*  init_filename = "fcc.xyz";
const int N_test = 500 ;

double density;
double beta; 


/* Simulation variables */
int n_particles = 0;
double e_cut;
double (*r)[NDIM];
double *size;
double box[NDIM];

double energy = 0.0;
double virial = 0.0;

double dummy;


/* Simulation variables *#include <sys/stat.h>
#include <sys/types.h>
/
int n_particles = 0;
double e_cut;
double (*r)[NDIM];
double *size;
double box[NDIM];

double energy = 0.0;
double virial = 0.0;

double dummy;



/*
#define N_lj 100 

double r_lj[N_lj];
double U[N_lj];
double F[N_lj];

double sigma = 1.0;
double epsilon = 1.0; 


void lenard_jones() {
    double base_c = sigma / r_cut;
    double e_cut = 4 * epsilon * (pow(base_c, 12) - pow(base_c, 6));

    for (int n = 0; n < N_lj; n++) {
        // Linear spacing from a small epsilon (to avoid r=0) up to r_cut
        // Starting at (n+1) prevents division by zero
        r_lj[n] = ((double)(n + 1) / N_lj) * r_cut;

        double r = r_lj[n];
        double s_over_r = sigma / r;
        double sr6 = pow(s_over_r, 6);
        double sr12 = sr6 * sr6;

        // Potential Energy (Shifted)
        U[n] = 4 * epsilon * (sr12 - sr6) - e_cut;

        // Force (Analytical Derivative)
        // F = 24 * eps * [2*(sig/r)^12 - (sig/r)^6] / r
        F[n] = (24.0 * epsilon / r) * (2.0 * sr12 - sr6);
    }
}
*/

double epsilon = 1.0; 

double calculate_force_times_dist(double r, int n) {
    if (r >= r_cut*r_cut) return 0.0; 

    double s_over_r = 1 / r;
    double sr6 = pow(s_over_r, 3);
    double sr12 = sr6 * sr6;
    
    return (48.0* epsilon) * ( sr12 - 0.5*sr6);
}


double calculate_potential(double r, int n) {

    if (r >= r_cut*r_cut) return 0.0; 

    double s_over_r = 1 / r;
    double sr6 = pow(s_over_r, 3);
    double sr12 = sr6 * sr6;
    
    return (4*epsilon) * ( sr12 - sr6) - e_cut;
}

typedef struct{
    double average_pressure;
    double mu_excess;
}measurement_t;

/* Functions */
measurement_t measure(void){
    measurement_t result;

    double force = 0;
    for(int h=0; h<n_particles; h++){


        for(int i=0; i < h; i++){

            double distance2 = 0;

            for(int j=0; j<3; j++){

                double dist = (r[h][j] - r[i][j]);

                if(dist>0.5*box[j]){
                    dist -= box[j];
                }
                else if(dist<-0.5*box[j]){
                    dist += box[j];
                }
                
                distance2 += dist*dist;
                    
            }
            
            force += calculate_force_times_dist(distance2, i); // THIS CAN BE MORE EFFICIENT DONT SQRT

            // double s_over_r = size[i] / distance;
            // double sr6 = pow(s_over_r, 6);
            // double sr12 = sr6 * sr6;
            

            // if (distance2 >= r_cut) force +=0;
            // else force +=(48.0* epsilon / distance2) * ( sr12 - 0.5*sr6);
        }
    }
    
    
    // dsfmt_genrand()

    double tot=0;

    for(int i=0; i < N_test;i++){
        double energys=0;
        
        double p[3];
        for(int k=0; k<3; k++){
            p[k] = dsfmt_genrand()*box[k];
        }

        for(int j=0;j<n_particles;j++){

            double distance2 = 0;

            for(int k=0; k<3; k++){

                double dist = (p[k] - r[j][k]);

                if(dist>0.5*box[k]){
                    dist -= box[k];
                }

                else if(dist<-0.5*box[k]){
                    dist += box[k];
                }
                
                distance2 += dist*dist;
                    
            }
        energys+= calculate_potential(distance2, j);
        }
        tot +=exp(-beta* energys);
    }

    tot/=N_test;

  

    double volume_box = pow(box[0],NDIM);

    result.average_pressure = (n_particles/volume_box)/ beta     + (force)/3/volume_box;

    result.mu_excess = - 1/beta * log(tot); 

    return result;
}


/*
void read_data(void){
    FILE* fp = fopen(init_filename, "r");
    int n, d;
    double dmin,dmax;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
        double diameter;
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
}
*/

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




typedef struct{
    double energy;
    double virial;
}particle_info_t;


particle_info_t particle_energy_and_virial(int pid){
    particle_info_t info;
    info.energy = 0.0;
    info.virial = 0.0;
    int n, d;
    for(n = 0; n < n_particles; ++n){
        if(n == pid) continue;
        double dist2 = 0.0;
        for(d = 0; d < NDIM; ++d){
            double min_d = r[pid][d] - r[n][d];
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];//?????
            dist2 += min_d * min_d;
        }

        if(dist2 <= r_cut * r_cut){
            double temp = 1.0 / (dist2 * dist2 * dist2);
            info.energy += 4.0 * temp * (temp - 1.0) - e_cut;
            info.virial += 24.0 * temp * (2.0 * temp - 1.0);
        }
    }

    return info;
}


int move_particle(void){
    int rpid = n_particles * dsfmt_genrand();


    particle_info_t info = particle_energy_and_virial(rpid);

    double old_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    particle_info_t new_info = particle_energy_and_virial(rpid);

    double dE = new_info.energy - info.energy;
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        virial += new_info.virial - info.virial;
        return 1;
    }

    for(d = 0; d < NDIM; ++d) r[rpid][d] = old_pos[d];

    return 0;
}

void write_data(int step){
    char buffer[128];

    sprintf(buffer, "./data/beta_%lf___T_%lf/coords_step%07d.dat",  n_particles/ pow(box[0],3) , 1/(beta*epsilon), step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", 1.0);
    }
    fclose(fp);
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

int main(int argc, char* argv[]){
    double base_c = 1*1 / r_cut*r_cut;
    double e_cut = 4 * epsilon * (pow(base_c, 6) - pow(base_c, 3));


    double betas[] = {0.5, 1.0, 2.0};
    double densitys[] = {0.1,0.3,0.5,0.6,0.7,0.8,0.9,1.0,1.1};

    
    // Calculate how many elements are in the array
    int big = sizeof(betas) / sizeof(betas[0]);
    int big2 = sizeof(densitys) / sizeof(densitys[0]);

    for (int i = 0; i < big; i++) {
        for(int j = 0; j < big2; j++){
            density= densitys[j];
            beta = betas[i];
            printf("Current beta: %.1f\n", beta);
            

            assert(delta > 0.0);

            e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));

            read_data();

            if(n_particles == 0){
                printf("Error: Number of particles, n_particles = 0.\n");
                return 0;
            }

            set_density();

            int d;
            for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box[d]);

            int step, n;
            for(n = 0; n < n_particles; ++n){
                particle_info_t info = particle_energy_and_virial(n);
                energy += info.energy;
                virial += info.virial;
            }
            energy *= 0.5;
            virial *= 0.5;

            size_t seed = time(NULL);
            dsfmt_seed(seed);

            double volume = 1.0;
            for(d = 0; d < NDIM; ++d) volume *= box[d];

            printf("Starting volume: %f\n", volume);
            printf("Starting energy: %f\n", energy);
            printf("Starting virial: %f\n", virial);
            printf("Starting seed: %lu\n", seed);

            char name[256];
            sprintf(name, "./data/beta_%lf___T_%lf", n_particles/ pow(box[0],3) , 1/(beta*epsilon));

            int result = mkdir(name ,0755);
            char save[128];

            sprintf(save, "./data/beta_%lf___T_%lf/measurements.dat",  n_particles/ pow(box[0],3) , 1/(beta*epsilon));
            FILE* fp = fopen(save, "w");


            int accepted = 0;
            for(step = 0; step < mc_steps; ++step){
                for(n = 0; n < n_particles; ++n){
                    accepted += move_particle();
                }

                measurement_t ms = measure();

                fprintf(fp, "%d\t%f\t%f\n", step, ms.average_pressure, ms.mu_excess);

                if(step % output_steps == 0){
                    printf("Step %d. Move acceptance: %f.\n",
                        step, (double)accepted / (n_particles * output_steps)
                    );
                    accepted = 0;
                    write_data(step);
                }
            }

            fclose(fp);
        }
    }

    return 0;
}

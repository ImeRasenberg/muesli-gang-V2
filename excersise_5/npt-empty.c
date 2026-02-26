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

/* Initialization variables */
const int mc_steps = 40000;
const int output_steps = mc_steps/200;
// wanting to save the converged state and the acceptance of both
double inf[40000/200][6];


const double packing_fraction = 0.55;
const double diameter = 1.0;    
const char* init_filename = "fcc.xyz";

// constants
const double betaP = 15;

double delta = 0.2;

double dV_m = 0.7;


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

// randome volume change
double dV;


double dummy;

int converged_move = 0;
int converged_vol = 0;




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




int change_volume(){

    dV = (dsfmt_genrand()-0.5)*2*dV_m;

    double V = box[0]*box[0]*box[0];
    
    double mult_fac = cbrt(V+dV)/box[0];

    double V_new=1;

    for(int i=0;i<3;i++){
        V_new *=box[i]*mult_fac;
    }
    
    double acc = fmin(1, exp(-betaP*dV + n_particles*log(V_new/V)));
    if (dsfmt_genrand()>acc){
        return 0;
    }


    double r_c[n_particles][3];

    for(int i=0; i<n_particles; i++){

        for(int j=0; j<3; j++){

            r_c[i][j] = r[i][j] *mult_fac; 
        }
    }

    for(int h=0; h<n_particles; h++){


        for(int i=0; i < h; i++){

            double distance = 0;

            for(int j=0; j<3; j++){

                double dist = (r_c[h][j] - r_c[i][j]);

                if(dist>0.5*box[j]*mult_fac){
                    dist -= box[j]*mult_fac;
                }
                else if(dist<-0.5*box[j]*mult_fac){
                    dist += box[j]*mult_fac;
                }
                
                distance += dist*dist;
                    
            }

            
            if (distance<(0.5*(size[h]+size[i]))*(0.5*(size[h]+size[i]))){
                // printf("volume cannot change there is overlap\n");
                return 0;
            }
        }

    }




    for(int i=0;i<3;i++){
        box[i]*=mult_fac;
    }
    for(int i=0; i<n_particles; i++){
        for(int j=0; j<3; j++){
            r[i][j] = r_c[i][j];
        }
    }

    // printf("volume changed\n");
    return 1;

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

    char name[256];
    sprintf(name, "./data/data %i", (int)floor(betaP));
    int result = mkdir(name ,0755);

    sprintf(buffer, "./data/data %i/coords_step%07d.dat",(int)floor(betaP), step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0, box[d]);
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


void info_2_file(){
    char new_name[128];
    sprintf(new_name, "./data/data %i/info.dat",(int)floor(betaP));

    FILE *print_coords; // inititialises a file variable
// char *new_name = "info.dat";
    print_coords = fopen(new_name,"w");

    double size_average = 0;
    for(int i = 0; i<n_particles;i++){
        size_average += size[i];
    }
    size_average/=n_particles;

    fprintf(print_coords, "%lf\t%lf\t%lf\t%i\t%lf\n", delta, dV_m, betaP, n_particles, size_average);

    for(int i=0;i<output_steps;i++){

        fprintf(print_coords, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", inf[i][0], inf[i][1],inf[i][2],inf[i][3],inf[i][4],inf[i][5]);
    }

    fclose(print_coords);


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
    int ind =0;
    int accepted_dv = 0;

    for(step = 1; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }
        accepted_dv += change_volume();

        if(step % output_steps == 0){

            double acceptance_move =  (double)accepted / (n_particles * output_steps);
            double acceptance_vol = (double)accepted_dv / (output_steps);

            printf("Step %d. Move acceptance: %lf.\n", step,  acceptance_move);
 
            printf("Step %d. Volume change acceptance: %lf.\n", step, acceptance_vol);

            write_data(step);

            if(converged_vol<4){
                if (acceptance_vol>0.55){
                    dV_m *= 1.1;
                }
                else if (acceptance_vol<0.45){
                    dV_m *= 0.9;
                }
                else{
                    converged_move =0;
                    converged_vol ++;
                }
            }

            if(converged_move<4){
                if (acceptance_move>0.55){
                    delta *= 1.1;
                }
                else if (acceptance_move<0.45){
                    delta *= 0.9;
                }
                else{
                    converged_move ++;

                }
            }

            inf[ind][2]=(double)acceptance_vol;
            inf[ind][4]=(double)acceptance_move;
            inf[ind][3]=(double)converged_vol;
            inf[ind][5]=(double)converged_move;
            inf[ind][0]=(double)step;
            inf[ind][1]=(double)box[0]*box[0]*box[0];
            ind++;


            

            accepted = 0;
            accepted_dv = 0;
        }

        
        
    }



    info_2_file();

    printf("done");
    return 0;
}

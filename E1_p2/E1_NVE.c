#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>
#include "../downloads/mt19937.h"

#include <sys/stat.h>
#include <sys/types.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512

const char*  init_filename = "fcc.xyz";
int n_particles = 0;
double (*r)[NDIM];
double (*F)[NDIM];
double (*Fn)[NDIM];
double (*v)[NDIM];
double *size;
double box[NDIM];
double dummy;

double dt = 5E-4; // the size of the timesteps
#define M 1000
double time_array[M];
double E[M];

double epsilon = 2/3;
double sigma = 1.0;
double r_cut = 2.5;
double e_cut = 0.0;
double density = 0.5;
double beta = 1.0;
double mass = 1.0;
double std;


double calculate_force_times_dist(double r) {
    if (r >= r_cut*r_cut) return 0.0; 

    double s_over_r = sigma*sigma / r;
    double sr6 = pow(s_over_r, 3);
    double sr12 = sr6 * sr6;
    
    return (48.0* epsilon) * ( sr12 - 0.5*sr6);
}
double calculate_potential(double r) {

    if (r >= r_cut*r_cut) return 0.0; 

    double s_over_r = sigma*sigma / r;
    double sr6 = pow(s_over_r, 3);
    double sr12 = sr6 * sr6;
    
    return (4*epsilon) * ( sr12 - sr6) - e_cut;
}
double energy;
void calc_forces(){
    for(int i=0; i<n_particles; i++) 
        for(int k=0; k<NDIM; k++) Fn[i][k] = 0.0;

    energy=0;

    double r_cut2 = r_cut*r_cut;
    for(int i = 0; i<n_particles; i++){
        for(int j = 0; j<i; j++){
            double r_ij[NDIM];
            double r2 = 0;
            
            for(int k=0; k < NDIM; k++) {
                r_ij[k] = r[i][k] - r[j][k];
                // printf("distance is %lf\n",r[i][k]);
                // Periodic Boundary Conditions
                if (r_ij[k] >  0.5 * box[k]) r_ij[k] -= box[k];
                if (r_ij[k] < -0.5 * box[k]) r_ij[k] += box[k];
                double len = r_ij[k] * r_ij[k];
                r2 += len;

                if(len < r_cut2){
                    double force = calculate_force_times_dist(len);
                    Fn[i][k] += force;
                    Fn[j][k] -= force;
                }
                
            }
            double energy_temp = calculate_potential(r2);
            // printf("distance is %lf\n",r2);/
            energy+= energy_temp; 
        }
    }

}

typedef struct{
    double N1[NDIM];
    double N2[NDIM];
}ran;
ran get_gaussian_nums(){
    ran ran_nums;
    for(int k=0; k<NDIM;k++){
        double U = dsfmt_genrand();
        // printf("U %lf\n",U);
        double V = dsfmt_genrand();
        // printf("V %lf\n",U);
        ran_nums.N1[k] = sqrt(-2.0*log(U))*cos(2*M_PI*V)*std;
        
        ran_nums.N2[k] = sqrt(-2.0*log(U))*sin(2*M_PI*V)*std;
    }
    return ran_nums;
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

void read_data(void){
   
    FILE *read_cords;
    read_cords = fopen(init_filename, "r");

    fscanf(read_cords, "%d\n", &n_particles); 

    r = malloc(n_particles * sizeof * r); 
    v = malloc(n_particles * sizeof * v);
    F = malloc(n_particles * sizeof * F);

    for(int i=0; i<n_particles; i++) 
        for(int k=0; k<NDIM; k++) F[i][k] = 0.0;

    Fn = malloc(n_particles * sizeof * Fn); 
    size = malloc(n_particles * sizeof * size); 

    for(int i = 0; i<NDIM; i++){
        fscanf(read_cords, "%lf %lf", &dummy, &box[i]);

    }

    for(int i = 0; i<n_particles; i++){
        fscanf(read_cords, "%lf %lf %lf %lf", &r[i][0], &r[i][1], &r[i][2], &size[i]);
    }

    fclose(read_cords);
}

int main(){
    std = 1/beta/mass;
    // printf("sdt %lf/n",std);
    e_cut = 4.0 * (pow(sigma / r_cut, 12.0) - pow(sigma / r_cut, 6.0));
    printf("sdt %lf/n",std);
    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_density();


    double dts[] = {1E-2,5E-2,1E-3,5E-3,1E-4,1E-5,1E-6,1E-7,5E-5,5E-6,5E-7};
    int big = sizeof(dts) / sizeof(dts[0]);

    for(int j=0;j<big;j++){
        dt = dts[j];
        size_t seed = time(NULL);
        dsfmt_seed(seed);


        for(int i=0;i< (int)floor(n_particles/2);i++){
            ran rannums = get_gaussian_nums();
            for(int k=0;k<NDIM;k++){
                v[(int)(2*i)][k] = rannums.N1[k];
                v[(int)(2*i)+1][k] = rannums.N2[k];
                // printf("speed is %lf\n ",v[(int)(2*i)][k]);
                // printf("speed is %lf\n ",v[(int)(2*i)+1][k]);
            }
        }
        for(int i=0; i<n_particles; i++) 
            for(int k=0; k<NDIM; k++) F[i][k] = 0.0;

        for(int c =0; c<M; c++){
            time_array[c] = (double)c*dt;
            
            for(int i=0;i<n_particles;i++){
                for(int j=0;j<NDIM;j++){
                    r[i][j] += v[i][j]*dt + Fn[i][j]/2/mass * dt*dt;
                    v[i][j] += F[i][j]/2/mass *dt;
                }
            }
            calc_forces();
            double E_kin = 0;

            for(int i=0;i<n_particles;i++){
                for(int j=0;j<NDIM;j++){
                    v[i][j] += Fn[i][j]/2/mass *dt;
                    F[i][j] = Fn[i][j];
                    
                    E_kin += 0.5*mass*v[i][j]*v[i][j];
                    // printf("speed squared is %lf\n",v[i][j]*v[i][j]);
                    // printf("mass is %lf\n",mass);
                    // printf("energy is %lf\n",E_kin);
                }
            }

            E[c] = energy+E_kin;
        }
        char filename[100];
        sprintf(filename, "data/energy_vs_time_dt_%.8f.txt", dt);

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# Time    Energy\n");
        for (int i = 0; i < M; i++) {
            fprintf(fp, "%f    %f\n", time_array[i], E[i]);
        }
        fclose(fp);
    }

}

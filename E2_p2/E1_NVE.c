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
#define N 256

const char*  init_filename = "fcc.xyz";
int n_particles = 0;
double (*r)[NDIM];
double (*F)[NDIM];
double (*v)[NDIM];
double *size;
double box[NDIM];
double dummy;

double dt = 1E-4; // the size of the timesteps
#define M 20000
double time_array[M];
double E[M][2];

double epsilon = 1;
double sigma = 1.0;
double r_cut = 0.0;
double e_cut = 0.0;
double density = 0.6;
double beta = 1.0;
double mass = 1.0;
double std;

double MSD; 
double MSD_arr[M];
double r_d[N][NDIM]; // displacement

double calculate_force_over_r(double r2) {
    if (r2 >= r_cut * r_cut) return 0.0; 

    double inv_r2 = 1.0 / r2;                 
    double sr2 = (sigma * sigma) * inv_r2;    
    double sr6 = sr2 * sr2 * sr2;             
    double sr12 = sr6 * sr6;                  
    
    // This returns F(r)/r
    return (48.0 * epsilon * inv_r2) * (sr12 - 0.5 * sr6);
}


double calculate_potential(double r2) {
    if (r2 >= r_cut*r_cut) return 0.0;

    double inv_r2 = 1.0 / r2;
    double sr2 = sigma*sigma * inv_r2;
    double sr6 = sr2*sr2*sr2;
    double sr12 = sr6*sr6;

    return 4*epsilon*(sr12 - sr6) - e_cut;
}


double energy;
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
    size = malloc(n_particles * sizeof * size); 

    for(int i = 0; i<NDIM; i++){
        fscanf(read_cords, "%lf %lf", &dummy, &box[i]);

    }

    for(int i = 0; i<n_particles; i++){
        fscanf(read_cords, "%lf %lf %lf %lf", &r[i][0], &r[i][1], &r[i][2], &size[i]);
    }

    fclose(read_cords);
}

void MeanSquaredDis (void){
    MSD = 0;
    for (int n = 0; n<n_particles; n++){
        double dis2 = 0;
        for (int d = 0; d<NDIM ; d++){
            dis2 += r_d[n][d]*r_d[n][d] ;
        }
    MSD += dis2;
    }
    MSD /= n_particles;
}

void write_MSD(void){
    char filename[100];
        sprintf(filename, "data/MD_MSD.txt");

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# Time    MSD\n");
        for (int t = 0; t < M; t++) {
            fprintf(fp, "%lf\t%lf\n", time_array[t], MSD_arr[t]);
        }
        fclose(fp);
}

int main(){
    printf("starting");
    r_cut = pow(2, 1.0/6.0)*sigma;
    std = sqrt(1/beta/mass);
    // printf("sdt %lf/n",std);
    e_cut = 4.0 * (pow(sigma / r_cut, 12.0) - pow(sigma / r_cut, 6.0));



    double dts[] = {1E-4};
    int big = sizeof(dts) / sizeof(dts[0]);

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    for(int o=0;o<big;o++){
        dt = dts[o];

        read_data();

        if(n_particles == 0){
            printf("Error: Number of particles, n_particles = 0.\n");
            return 0;
        }

        set_density();

        for(int i=0;i< (int)floor(n_particles/2);i++){
            ran rannums = get_gaussian_nums();
            for(int k=0;k<NDIM;k++){
                v[(int)(2*i)][k] = rannums.N1[k];
                v[(int)(2*i)+1][k] = rannums.N2[k];
                // printf("speed is %lf\n ",v[(int)(2*i)][k]);
                // printf("speed is %lf\n ",v[(int)(2*i)+1][k]);
            }
        }

        for(int k=0;k<NDIM;k++){
        double v_cm = 0;
        for(int i=0;i<n_particles;i++) v_cm += v[i][k];
        v_cm /= n_particles;

        for(int i=0;i<n_particles;i++) v[i][k] -= v_cm;
        }

        calc_forces();

        for(int c = 0; c < M; c++) {
            time_array[c] = (double)c * dt;


            for(int i = 0; i < n_particles; i++) {
                for(int j = 0; j < NDIM; j++) {
                    double dx = v[i][j] * dt + 0.5 /mass * F[i][j] * dt *dt;
                    r[i][j] += dx;
                    r_d[i][j]+= dx;
                    if (r[i][j] < 0) r[i][j] += box[j];
                    if (r[i][j] >= box[j]) r[i][j] -= box[j];

                    v[i][j] += (F[i][j] / (2.0 * mass)) * dt;
                }
            }


            calc_forces(); 


            double E_kin = 0;
            for(int i = 0; i < n_particles; i++) {
                double v2_sum = 0;
                for(int j = 0; j < NDIM; j++) {

                    v[i][j] += (F[i][j] / (2.0 * mass)) * dt;
                    v2_sum += v[i][j] * v[i][j];
                }
                E_kin += 0.5 * mass * v2_sum;
            }

            E[c][0] = energy; 
            E[c][1] = E_kin;
            MeanSquaredDis();
            MSD_arr[c]= MSD;
        }
        write_MSD();
        char filename[100];
        sprintf(filename, "data/MD_energy.txt");

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# Time    Energy\n");
        for (int i = 0; i < M; i++) {
            fprintf(fp, "%lf\t%lf\t%lf\n", time_array[i], E[i][0],E[i][1]);
        }
        fclose(fp);
    }
    printf("done");

}

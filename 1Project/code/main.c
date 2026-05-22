//  cd /home/ime_rasenberg/Documents/muesli-gang-V2/1Project/code/
//  gcc main.c -o main -lm
//  ./main

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mt19937.h"

const double pi = 3.14159265358979323846;

#define N 20
#define M 30000

double spin[N][N][5]; // The lattice in which all spins reside with their angles

// constants of the Hamiltonian
double J = 1;
double D = 1;
double Hz = 0.5;

// Constants of the simmulation
double beta = 0.1;
double dang = 0.5*pi;

// Simmulation metrics
int accepted = 0;
double Energy=0;



void generate_random_spin(){

    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){

            double phi = pi * dsfmt_genrand();
            double psi = 2.0 * pi * dsfmt_genrand();


            spin[i][j][0] = sin(phi) * cos(psi);
            spin[i][j][1] = sin(phi) * sin(psi);
            spin[i][j][2] = cos(phi);
            spin[i][j][3] = phi;
            spin[i][j][4] = psi;

        }
    }
}

void WriteState2File(){
    char filename[100];
    sprintf(filename, "Data/Spin_orientiation.txt");

    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("ERROR: Could not create file! Make sure a folder named 'Data' exists in this directory.\n");
        return; 
    }

    fprintf(fp, "%d\n", N);

    for(int k = 0; k < 3; k++){
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                fprintf(fp, "%lf\t", spin[i][j][k]);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

double dot_prod(double A[5], double B[5]){
    double res=0;
    for(int k =0; k<3; k++){
        res+=A[k]*B[k];
    }
    return res;
}

double cross_x(double A[5], double B[5]){
    return A[1]*B[2] - A[2]*B[1]; 
}
double cross_y(double A[5], double B[5]){
    return A[2]*B[0] - A[0]*B[2]; 
}


double get_energy(double spin[N][N][5]){
    double Energy_n = 0;
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            // interaction choice
            int i_2 = i;
            if (i == N) i_2=0;
            int j_2 = j;
            if (j == N) j_2=0;

            // Magnetic Field energy
            double E_H = Hz * spin[i][j][2];

            // spin spin alinging interaction
            double E_J = 0;
            E_J += -J*dot_prod(spin[i][j],spin[i_2][j]);

            E_J += -J*dot_prod(spin[i][j],spin[i][j_2]);

            // spin orbit coupling
            double E_D = 0;
            E_D += -D*cross_x(spin[i][j],spin[i][j_2]);
            E_D += -D*cross_y(spin[i][j],spin[i][j_2]);

            E_D += -D*cross_x(spin[i][j],spin[i_2][j]);
            E_D += -D*cross_y(spin[i][j],spin[i_2][j]);

            Energy_n += E_H + E_J + E_D;
        }
    }

    return Energy_n;
}

int change_particle(){
    double spin_n[N][N][5];
    int n_1 = floor(N*dsfmt_genrand());
    int n_2 = floor(N*dsfmt_genrand());
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            for(int k=0; k<5; k++){
                spin_n[i][j][k] = spin[i][j][k];
            }
        }
    }

    // bounded
    spin_n[n_1][n_2][3] += dang * dsfmt_genrand(); // domain 0,pi
    if  (spin_n[n_1][n_2][3] > pi) spin_n[n_1][n_2][3] -= pi;
    else if (spin_n[n_1][n_2][3] < 0) spin_n[n_1][n_2][3] += pi;
    spin_n[n_1][n_2][4] += dang * dsfmt_genrand(); // domain 0, 2 pi

    spin_n[n_1][n_2][0] = sin(spin_n[n_1][n_2][3]) * cos(spin_n[n_1][n_2][4]);
    spin_n[n_1][n_2][1] = sin(spin_n[n_1][n_2][3]) * sin(spin_n[n_1][n_2][4]);
    spin_n[n_1][n_2][2] = cos(spin_n[n_1][n_2][3]);


    // double phi = pi * dsfmt_genrand();
    // double psi = 2.0 * pi * dsfmt_genrand();


    // spin[n_1][n_2][0] = sin(phi) * cos(psi);
    // spin[n_1][n_2][1] = sin(phi) * sin(psi);
    // spin[n_1][n_2][2] = cos(phi);


    double dE = get_energy(spin_n) - Energy;

    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        return 0;
    }

    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            for(int k=0; k<3; k++){
                spin[i][j][k] = spin_n[i][j][k];
            }
        }
    }
    Energy += dE;

    return 1; 

}

int main(void){
    dsfmt_seed(time(NULL));

    generate_random_spin();
    WriteState2File();

    Energy = get_energy(spin);
    // printf("initial energy is %lf\n", Energy);


    char save[128];
    sprintf(save, "Data/Energy_steps.txt");
    FILE* fp = fopen(save, "w");


    for(int count=0; count<M; count++){
        accepted += change_particle();
        fprintf(fp, "%d\t%lf\t%lf\n", count, Energy, accepted/(double)count);
    }
    fclose(fp);

    WriteState2File();

    return 0;
}

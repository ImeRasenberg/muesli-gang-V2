//  cd /home/ime_rasenberg/Documents/muesli-gang-V2/1Project/code/
//  gcc main.c -o main -lm
//  ./main

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mt19937.h"

const double pi = 3.14159265358979323846;

#define N 20
#define M 1e5

double spin[N][N][3]; // The lattice in which all spins reside with their angles
double spin_n[N][N][3]; // The updated lattice

// constants of the Hamiltonian
double J = 1;
double D = 1;
double Hz = 0.5;

// Constants of the simmulation
double beta = 1;
double dang = 1;
int NC = 1;
int n1 = 0;
int n2 = 0;

// Simmulation metrics
int accepted = 0;
double Energy=0;



void generate_random_spin(){

    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){

            double z   = 2.0 * dsfmt_genrand() - 1.0;
            double psi = 2.0 * pi * dsfmt_genrand();

            double r = sqrt(1.0 - z*z);

            spin[i][j][0] = r * cos(psi);
            spin[i][j][1] = r * sin(psi);
            spin[i][j][2] = z;
        }
    }
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            for(int k=0; k<3; k++){
                spin_n[i][j][k] = spin[i][j][k];
            }
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

double dot_prod(double A[3], double B[3]){
    double res=0;
    for(int k =0; k<3; k++){
        res+=A[k]*B[k];
    }
    return res;
}
double cross_x(double A[3], double B[3]){
    return A[1]*B[2] - A[2]*B[1]; 
}
double cross_y(double A[3], double B[3]){
    return A[2]*B[0] - A[0]*B[2]; 
}


double get_energy_tot(double spin[N][N][3]){
    double Energy_n = 0;
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            // interaction choice
            int i_2 = (i + 1) % N;
            int j_2 = (j + 1) % N;

            // Magnetic Field energy
            double E_H = -Hz * spin[i][j][2];

            // spin spin alinging interaction
            double E_J = 0;

            E_J += -J*dot_prod(spin[i][j],spin[i_2][j]);
            E_J += -J*dot_prod(spin[i][j],spin[i][j_2]);

            // spin orbit coupling
            double E_D = 0;
            E_D += -D*cross_x(spin[i][j],spin[i_2][j]);
            E_D += -D*cross_y(spin[i][j],spin[i][j_2]);

            Energy_n += E_H + E_J + E_D;
        }
    }

    return Energy_n;
}


// THIS IS A QUESTION FOR THE TA'S, WHY ISNT THIS EXACTLY DELTA_E???
double get_energy_n(double spin[N][N][3]){
    double Energy_n = 0;

    for (int di = -NC; di <= NC; di++) {
        int i = (n1 + di + N) % N;
        for (int dj = -NC; dj <= NC; dj++) {
            int j = (n2 + dj + N) % N;

            // Magnetic Field energy
            double E_H = -Hz * spin[i][j][2];

            // indexes for spin spin interactions
            int i_2 = (i + 1 + N) % N;
            int i_3 = (i - 1 + N) % N;
            int j_2 = (j + 1 + N) % N;
            int j_3 = (j - 1 + N) % N;


            // spin spin alinging interaction
            double E_J = 0;
            E_J += -J/2*dot_prod(spin[i][j],spin[i_2][j]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i_3][j]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i][j_2]);
            E_J += -J/2*dot_prod(spin[i][j],spin[i][j_3]);

            // spin orbit coupling
            double E_D = 0;
            E_D += -D/2*cross_x(spin[i][j],spin[i_2][j]);
            E_D += -D/2*cross_x(spin[i][j],spin[i_3][j]);
            E_D += -D/2*cross_y(spin[i][j],spin[i][j_2]);
            E_D += -D/2*cross_y(spin[i][j],spin[i][j_3]);

            Energy_n += E_H + E_J + E_D;
        }
    }

    return Energy_n;
}

int change_particle(){
    spin_n[n1][n2][0] += dang*(2.0*dsfmt_genrand()-1.0);
    spin_n[n1][n2][1] += dang*(2.0*dsfmt_genrand()-1.0);
    spin_n[n1][n2][2] += dang*(2.0*dsfmt_genrand()-1.0);

    double norm = sqrt(
        spin_n[n1][n2][0]*spin_n[n1][n2][0] +
        spin_n[n1][n2][1]*spin_n[n1][n2][1] +
        spin_n[n1][n2][2]*spin_n[n1][n2][2]
    );

    spin_n[n1][n2][0] /= norm;
    spin_n[n1][n2][1] /= norm;
    spin_n[n1][n2][2] /= norm;


    double dE = get_energy_n(spin_n) - get_energy_n(spin);

    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){

        // for(int i=0; i < N; i++){
        //     for(int j=0; j < N; j++){
        //         for(int k=0; k<3; k++){
        //             spin[i][j][k] = spin_n[i][j][k];
        //         }
        //     }
        
        // }
        spin[n1][n2][0] = spin_n[n1][n2][0];
        spin[n1][n2][1] = spin_n[n1][n2][1];
        spin[n1][n2][2] = spin_n[n1][n2][2];
        Energy += dE;
        return 1;
    }

    spin_n[n1][n2][0] = spin[n1][n2][0];
    spin_n[n1][n2][1] = spin[n1][n2][1];
    spin_n[n1][n2][2] = spin[n1][n2][2];
    return 0; 

}

int main(void){
    dsfmt_seed(time(NULL));

    generate_random_spin();
    WriteState2File();

    Energy = get_energy_tot(spin);
    printf("initial energy is %lf\n", Energy);


    char save[128];
    sprintf(save, "Data/Energy_steps.txt");
    FILE* fp = fopen(save, "w");

    for(int count=1; count<M+1; count++){
        n1 = floor(N*dsfmt_genrand());
        n2 = floor(N*dsfmt_genrand());
        accepted += change_particle();
        fprintf(fp, "%d\t%lf\t%lf\n", count, Energy, accepted/(double)count);
    }
    fclose(fp);

    WriteState2File();

    return 0;
}

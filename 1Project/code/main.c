//  cd /home/ime_rasenberg/Documents/muesli-gang-V2/1Project/code/
//  gcc main.c -o main -lm
//  ./main

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mt19937.h"

const double pi = 3.14159265358979323846;

#define N 40
#define M 16e4

double spin[N][N][3]; // The lattice in which all spins reside with their angles
double spin_n[N][N][3]; // The updated lattice

// constants of the Hamiltonian
double J = 1;
double D = 0.2;
double Hz = 0.08;

// Constants of the simmulation
double beta = 5;
double dang = 1.0;
int NC = 1;
int n1 = 0;
int n2 = 0;

// Simmulation metrics
int accepted = 0;
double Energy=0;
double Q = 0;



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
double get_delta_energy(int i, int j) {
    // Calculate the change in the spin vector
    double dS[3];
    for(int k=0; k<3; k++) {
        dS[k] = spin_n[i][j][k] - spin[i][j][k];
    }

    // Neighbors with periodic boundary conditions
    int ip = (i + 1) % N;
    int im = (i - 1 + N) % N;
    int jp = (j + 1) % N;
    int jm = (j - 1 + N) % N;

    // Magnetic Field contribution
    double dE_H = -Hz * dS[2];

    // Heisenberg Exchange contribution (J)
    double dE_J = -J * (
        dot_prod(dS, spin[ip][j]) +
        dot_prod(dS, spin[im][j]) +
        dot_prod(dS, spin[i][jp]) +
        dot_prod(dS, spin[i][jm])
    );

    // Dzyaloshinskii-Moriya Interaction (DMI) contribution (D)
    double dE_D = 0;
    
    dE_D += -D * cross_x(dS, spin[ip][j]);          // Interaction where (i,j) is the left actor
    dE_D += -D * cross_x(spin[im][j], dS);          // Interaction where (i,j) is the right actor

    // Along Y direction: cross_y(S_j, S_j+1)
    dE_D += -D * cross_y(dS, spin[i][jp]);          // Interaction where (i,j) is the bottom actor
    dE_D += -D * cross_y(spin[i][jm], dS);          // Interaction where (i,j) is the top actor

    return dE_H + dE_J + dE_D;
}
int change_particle(){
    // changing the spin with some random ammount
    spin_n[n1][n2][0] += dang*(2.0*dsfmt_genrand()-1.0);
    spin_n[n1][n2][1] += dang*(2.0*dsfmt_genrand()-1.0);
    spin_n[n1][n2][2] += dang*(2.0*dsfmt_genrand()-1.0);

    double norm = sqrt(
        spin_n[n1][n2][0]*spin_n[n1][n2][0] +
        spin_n[n1][n2][1]*spin_n[n1][n2][1] +
        spin_n[n1][n2][2]*spin_n[n1][n2][2]
    );

    spin_n[n1][n2][0] /= norm;
    spin_n[n1][n2][1] /= norm; //completely solid white background everywhere else, showing that the rest of the lattice has perfectly aligned spins (qij​=0).norm;
    spin_n[n1][n2][2] /= norm;

    // The change in enery
    double dE = get_delta_energy(n1 , n2);

    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
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


void cross_prod(double A[3], double B[3], double C[3]){
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}
double triangle_charge(double S1[3], double S2[3], double S3[3]){
    double AxB[3];

    cross_prod(S2, S3, AxB);
    double num = S1[0]*AxB[0] + S1[1]*AxB[1] + S1[2]*AxB[2];

    double d12 = S1[0]*S2[0] + S1[1]*S2[1] + S1[2]*S2[2];
    double d23 = S2[0]*S3[0] + S2[1]*S3[1] + S2[2]*S3[2];
    double d31 = S3[0]*S1[0] + S3[1]*S1[1] + S3[2]*S1[2];

    double den = 1.0 + d12 + d23 + d31;

    return 2.0 * atan2(num, den) / (4.0 * pi);
}
double get_Q(){
    double Q = 0.0;
    double CD[N][N];

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){

            int ip = (i + 1) % N;
            int jp = (j + 1) % N;

            double S_ij[3]   = {spin[i][j][0],   spin[i][j][1],   spin[i][j][2]};
            double S_ipj[3]  = {spin[ip][j][0],  spin[ip][j][1],  spin[ip][j][2]};
            double S_ijp[3]  = {spin[i][jp][0],  spin[i][jp][1],  spin[i][jp][2]};
            double S_ipjp[3] = {spin[ip][jp][0], spin[ip][jp][1], spin[ip][jp][2]};

            CD[i][j] = 0;
            // two triangles per plaquette
            CD[i][j] += triangle_charge(S_ij, S_ipj, S_ijp);
            CD[i][j] += triangle_charge(S_ipj, S_ipjp, S_ijp);

            Q += CD[i][j];
        }
    }

    return Q;
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

    // we start hot and then cool down
    int tot =0;
    beta = 0.5;
    double betas[] = {0.5, 1.0, 2.0, 4.0};
    int big = sizeof(betas) / sizeof(betas[0]);

    for(int s = 0; s<big; s++){
        beta = betas[s];
        for(int count=1; count<M+1; count++){
            if(count % 1000 == 0){
                Q = get_Q();
            }
            n1 = floor(N*dsfmt_genrand());
            n2 = floor(N*dsfmt_genrand());
            accepted += change_particle();
            tot+=1;
            fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", tot, Energy, accepted/(double)tot, beta, Q);
        }
    }

    fclose(fp);
    WriteState2File();

    return 0;
}

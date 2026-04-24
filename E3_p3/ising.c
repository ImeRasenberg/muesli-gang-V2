#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 40
#define M 10000
double J = 1; //coupeling constant
double beta = 0.33;

double energy; // the energy of the system
double magnet;  // the magnetization of the system
double heat_cp; // the heat capacity of the system

int S[N][N]; // the spin array
double energy_arr[M];
double magnet_arr[M];
double heat_cp_arr[M];

void fill_S (void){
    /*setting everyting spin up*/
    for (int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            S[i][j] = 1;
        }
    }
}

void print_lattice(void) {
    /*function that prints the grid */
    printf("Current spin lattice:\n");
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            printf("%d ", S[r][c]);
        }
        printf("\n");
    }
}

double calculate_energy (void){
    /*calcuates and retuns the energy based on the current S grid*/
    double E = 0;
    for (int r = 0; r <N ; r++){ // sum over rows coloms
        for (int c = 0; c<N; c++){
        for (int dr = -1; dr<=1; dr++ ){ // sum over neighourbs
            for (int dc = -1; dc<=1 ; dc++){
                if (dr == 0 && dc == 0) continue;
                
                int r_n = r +dr; // neighbouring row
                int c_n = c +dc; // neighbouring collom

                if (r_n >= N) r_n -= N; // apply periodic boundary conditions
                if (r_n < 0) r_n += N;
                
                if (c_n >= N) c_n -= N;
                if (c_n < 0) c_n += N;

                E += -J/2.0 * S[r][c] * S[r_n][c_n];
             }
          }
        
        }
    }
    return E;
}


double delta_E(int r, int c) {
    double sum = 0;
    for (int dr = -1; dr <= 1; dr++) {
        for (int dc = -1; dc <= 1; dc++) {
            if (dr == 0 && dc == 0) continue;
            int rn = (r + dr + N) % N;
            int cn = (c + dc + N) % N;
            sum += S[rn][cn];
        }
    }
    
    return 2* J * S[r][c] * sum;
}

double calculate_magnetization(void){
     /*calcuates and retuns the magnetization based on the current S grid*/
    double Mag = 0;
    for (int r = 0; r<N; r++){
        for (int c = 0; c<N; c++){
            Mag += S[r][c];
        }
    }
    return Mag;
}

int Monte_carlo (void){

    int r = (int)(N * dsfmt_genrand()); // pick random spin
    int c = (int)(N * dsfmt_genrand());
  
    double dE= delta_E(r, c);
    double acc = exp(- beta * dE);

    if (dsfmt_genrand()< acc){ // accept trail move and flip spin, add energy change and magnetization change
        S[r][c] *= -1;
        energy += dE;
        magnet += 2* S[r][c];
        return 1; 
    }

    return 0; 

}

void write_data(void){
    char filename[100];
        sprintf(filename, "data/energy_magnet_vs_time%.5f.txt", beta);

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# step,  Energy, magenetization \n");
        for (int t = 0; t < M; t++) {
            fprintf(fp, "%d\t%lf\t%lf\n", t, energy_arr[t],magnet_arr[t]);
        }
        fclose(fp);
}


int main (void){
    dsfmt_seed(time(NULL));
    fill_S();

    //print_lattice();

    energy = calculate_energy();
    magnet = calculate_magnetization();

    energy_arr[0]= energy;
    magnet_arr[0]= magnet;

    printf("the start energy is %lf\n", energy);
    printf("the start magnetization is %lf\n", magnet);

    for (int step = 1; step <M; step ++){
        for (int i = 0; i< N* N; i++){
             Monte_carlo();
        }
        energy_arr[step]= energy;
        magnet_arr[step]= magnet;
       
    }
    
    double calc_energy = calculate_energy();
    double calc_mag = calculate_magnetization();

    printf("the final energy is %lf\n", energy);
    printf("the final calcuated energy is %lf\n", energy);

    printf("the final magnetization is %lf\n", magnet);
    printf("the final calcuated magnetization is %lf\n", calc_mag);

    //print_lattice();
    printf("done");
    write_data();
    return 0;
}




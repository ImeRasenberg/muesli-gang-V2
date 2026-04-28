#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 32
#define M 1000
#define MAX_T 14      // amount of temps (not clean but idk)
#define cutof 200 // cutof ater which we sample
#define INTERVAL (M - cutof)

double J = 1; //coupeling constant
double beta = 1.0;
double Temp = 1.0;


double energy; // the energy of the system
double magnet;  // the magnetization of the system
int S[N][N]; // the spin array
double magnet_arr[M];
double correlations[MAX_T][INTERVAL]; //the correlations




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
void calculate_correlation(double corr[][INTERVAL], int idx){
   
    int interval = M - cutof;
    int max_t = interval /2;

    double M_tavg = 0.0, M2_tavg = 0.0;
    
    for (int i = cutof; i< M ; i++){
        double m = magnet_arr[i];
        M_tavg += m;
        M2_tavg += m* m;
    }
    
    M_tavg /= interval;
    M2_tavg /= interval;

    double variance = M2_tavg - M_tavg * M_tavg ;

    for (int t = 0; t < max_t ; t++) {
       
        double chi = 0;
        int n_pairs = 0;

        for (int tau = cutof; tau < M- t; tau++){
            chi += magnet_arr[tau]*magnet_arr[t + tau];
            n_pairs++;
            }
        chi /= n_pairs;
        if (variance < 1e-8) {
        corr[idx][t] = (t == 0) ? 1.0 : 0.0;
            } 
        else {
            chi = (chi - M_tavg * M_tavg)/variance;
            corr[idx][t] = chi;
             }
        }
    }
void write_data(int idx, double Temp ){
    char filename[100];
        sprintf(filename, "data/correlation_vs_time_Temp:%.3f.txt", Temp);

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# time,  correlation,  \n");
        for (int i = 0; i < INTERVAL/2; i++) {
            fprintf(fp, "%d\t%lf\t%lf\n", i, correlations[idx][i], Temp);
        }
        fclose(fp);
}

int main (void){

    int interval = M- cutof;
    dsfmt_seed(time(NULL));
    clock_t end, start;
    double cpu_time_used;
    start = clock();
    
    fill_S();

    double T_list[] = {1.0,  1.5,  2.0, 2.5, 3.0, 3.5, 4, 4.5, 5, 5.3,  5.5, 6.0, 7.0, 8.0 }; // loop over Temps's
    int len_T = sizeof(T_list)/sizeof(T_list[0]);
 
    printf("temp runs %d\n", len_T);
    for (int idx = 0; idx < len_T; idx++) {
       
        Temp = T_list[idx];
        beta = 1.0/Temp;

        fill_S();
        energy = calculate_energy();
        magnet = calculate_magnetization();

        for (int step = 0; step <M; step ++){
            for (int i = 0; i< N* N; i++){
             Monte_carlo();
            }
             magnet_arr[step]= magnet/ (N*N);
        } 
        calculate_correlation(correlations, idx);
    }

    for (int i = 0; i<len_T; i++){
        Temp = T_list[i];
        write_data(i, Temp);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("done in %f seconds\n", cpu_time_used);
    
    return 0;
}




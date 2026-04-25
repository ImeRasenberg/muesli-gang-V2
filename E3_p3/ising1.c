#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 100
#define M 1000

int cutof = 100; // cutof ater which we sample
double J = 1; //coupeling constant
double beta = 1.0;
double Temp = 1.0;


double energy; // the energy of the system
double magnet;  // the magnetization of the system
double heat_cp; // the heat capacity of the system

int S[N][N]; // the spin array
double energy_arr[M];
double magnet_arr[M];
double heat_cp_arr[M];

double E;  // the averages
double E2;
double Cv;
double av_magnet;



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
        sprintf(filename, "data/energy_magnet_vs_time_beta:%.3f.txt", beta);

        FILE *fp = fopen(filename, "w");
        fprintf(fp, "# step,  Energy, magenetization \n");
        for (int t = 0; t < M; t++) {
            fprintf(fp, "%d\t%lf\t%lf\t%lf\n", t, energy_arr[t],magnet_arr[t], beta);
        }
        fclose(fp);
}

void calculate (void){
    E = 0;
    E2 = 0;
    Cv = 0;
    av_magnet =0;

    for (int step = cutof; step<M; step++){
        E += energy_arr[step];
        E2 += energy_arr[step]* energy_arr[step];
        av_magnet += magnet_arr[step];
    }

    E /= (M-cutof);
    E2 /= (M-cutof);
    av_magnet/= (M-cutof);
    Cv=  beta * beta *( E2 - E*E)/(N*N);

    E /= N*N;
    E2 /= (N*N)*(N*N);
    av_magnet/= N*N;
}



int main (void){

    dsfmt_seed(time(NULL));
    

    clock_t end, start;
    double cpu_time_used;
    start = clock();

    double T_list[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 , 1.8, 1.9, 2.0, 2.1, 2.2 , 2.3 , 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7 , 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8 ,4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7 , 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5 ,6.6, 6.7, 6.8, 6.9, 7.0 }; // loop over Temps's
    int len_T = sizeof(T_list)/sizeof(T_list[0]);
    double E_avg[len_T];
    double E2_avg[len_T];
    double Cv_arr[len_T];
    double av_magnet_arr[len_T];

   
  


    printf("temp runs %d\n", len_T);
    for (int idx = 0; idx < len_T; idx++) {
       
        Temp = T_list[idx];
        beta = 1/Temp;

        fill_S();
        energy = calculate_energy();
        magnet = calculate_magnetization();
        energy_arr[0]= energy;
        magnet_arr[0]= magnet;

        for (int step = 1; step <M; step ++){
            for (int i = 0; i< N* N; i++){
             Monte_carlo();
            }
        energy_arr[step]= energy;
        magnet_arr[step]= magnet;
        } 

        calculate();
        E_avg[idx] = E;
        E2_avg[idx]= E2;
        Cv_arr[idx] = Cv;
        av_magnet_arr[idx] = av_magnet;
  

        printf("temp, magnetization, energy, energy^2 , Cv %lf\t%lf\t%lf\t%lf\t%lf\n", Temp, av_magnet, E, E2, Cv);
        //write_data(); // write data for every beta
    }
    char filename[100];
    sprintf(filename, "data/averages.txt");
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "# temp,  magnetization, average energy, average squared energy, heat_capacity \n");
    for (int t = 0; t < len_T; t++) {
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", T_list[t],av_magnet_arr[t], E_avg[t], E2_avg[t], Cv_arr[t]);
        }
    fclose(fp);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("done in %f seconds\n", cpu_time_used);
    
    return 0;
}




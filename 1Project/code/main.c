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
double CD[N][N]; //the wrapping number in space
double CD_smooth[N][N]; //the wrapping number smoothed out

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


// starting the simmulation
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

// writing out the state into a file for nice plotting
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

// The functions that determain the MC simmulations
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

// Getting the Q-number and tracing places of skyrmions
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

// Tracing the skyrmions
void gaussian_filter_CD(double sigma){

    // Radius ≈ 3 sigma (same practical cutoff scipy uses)
    int radius = (int)round(4.0 * sigma);
    int size = 2 * radius + 1;

    double kernel[size][size];
    double norm = 0.0;

    // Build Gaussian kernel
    for(int dx = -radius; dx <= radius; dx++){
        for(int dy = -radius; dy <= radius; dy++){

            double r2 = dx*dx + dy*dy;

            kernel[dx + radius][dy + radius]
                = exp(-r2 / (2.0 * sigma * sigma));

            norm += kernel[dx + radius][dy + radius];
        }
    }

    // Normalize kernel
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            kernel[i][j] /= norm;
        }
    }

    // Convolution with periodic boundary conditions (mode="wrap")
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){

            double sum = 0.0;

            for(int dx = -radius; dx <= radius; dx++){
                for(int dy = -radius; dy <= radius; dy++){

                    int ii = (i + dx + N) % N;
                    int jj = (j + dy + N) % N;

                    sum += CD[ii][jj]
                         * kernel[dx + radius][dy + radius];
                }
            }

            CD_smooth[i][j] = sum;
        }
    }
}
#define MAX_PEAKS (N*N)
int max_count = 0;
int min_count = 0;
int max_x[MAX_PEAKS];
int max_y[MAX_PEAKS];
int min_x[MAX_PEAKS];
int min_y[MAX_PEAKS];
void desission_peak(int step, char f2[128]){

    max_count = 0;
    min_count = 0;

    int radius = 2;
    int Q_target = (int)round(Q);

    double largest_peak = 0.0;

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){

            double val = fabs(CD_smooth[i][j]);

            if(val > largest_peak){
                largest_peak = val;
            }
        }
    }

    double best_frac = 0.25;
    int Q_found = 0;

    for(int step = 0; step < 40; step++){

        // gradually relax threshold
        double frac = 0.25 - step * (0.23 / 39.0);

        if(frac < 0.02){
            frac = 0.02;
        }

        double threshold = frac * largest_peak;

        max_count = 0;
        min_count = 0;

        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){

                double center = CD_smooth[i][j];

                // threshold cut
                if(fabs(center) < threshold){
                    continue;
                }

                int is_max = 1;
                int is_min = 1;

                for(int dx = -radius; dx <= radius; dx++){

                    for(int dy = -radius; dy <= radius; dy++){

                        if(dx == 0 && dy == 0){
                            continue;
                        }

                        int ii = (i + dx + N) % N;
                        int jj = (j + dy + N) % N;

                        double neighbor = CD_smooth[ii][jj];

                        // maxima test
                        if(center < neighbor){
                            is_max = 0;
                        }

                        // minima test
                        if(center > neighbor){
                            is_min = 0;
                        }

                        // early exit
                        if(!is_max && !is_min){
                            goto done_check;
                        }
                    }
                }

done_check:

                if(is_max){
                    max_x[max_count] = i;
                    max_y[max_count] = j;
                    max_count++;
                }

                if(is_min){
                    min_x[min_count] = i;
                    min_y[min_count] = j;
                    min_count++;
                }
            }
        }

        Q_found = max_count - min_count;

        best_frac = frac;

        // success condition
        if(Q_found == Q_target){
            break;
        }
    }

    FILE *fp = fopen(f2, "a");
    if (fp == NULL) {
        return;
    }

    // 1. Start the main object and print basic simulation metrics
    fprintf(fp, "{\n");
    fprintf(fp, "  \"step\": %d,\n", step);
    fprintf(fp, "  \"Q\": %d,\n", Q_found);
    fprintf(fp, "  \"threshold_fraction\": %.6f,\n", best_frac);
    fprintf(fp, "  \"E\": %.6f,\n", Energy);
    fprintf(fp, "  \"J\": %.6f,\n", J);
    fprintf(fp, "  \"D\": %.6f,\n", D);
    fprintf(fp, "  \"Hz\": %.6f,\n", Hz);
    fprintf(fp, "  \"beta\": %.6f,\n", beta);
    

    // 2. Print positive peaks ("N+")
    fprintf(fp, "  \"N+\": {");
    for (int k = 0; k < max_count; k++) {
        fprintf(fp, "\"%d\": [%d, %d]", k, max_x[k], max_y[k]);
        if (k < max_count - 1) {
            fprintf(fp, ", "); // Comma between peak objects
        }
    }
    fprintf(fp, "},\n"); // Close N+ and add a comma for the next key

    // 3. Print negative peaks ("N-")
    fprintf(fp, "  \"N-\": {");
    for (int k = 0; k < min_count; k++) {
        fprintf(fp, "\"%d\": [%d, %d]", k, min_x[k], min_y[k]);
        if (k < min_count - 1) {
            fprintf(fp, ", "); // Comma between peak objects
        }
    }
    fprintf(fp, "}\n"); // Close N- (no trailing comma here because it's the last item!)

    // 4. Close the main object
    fprintf(fp, "}\n");

    fclose(fp);
}



// where the fun happends
int main(void){
    dsfmt_seed(time(NULL));

    generate_random_spin();
    WriteState2File();

    // stabalisation constants
    int window_size = N*N * 3; 
    double tolerance = 50;
    int number_windows = 2; // number of requiered stable windows



    int O = 30;
    for(int count1 = 0; count1<O+1; count1++){
        D=2/(double)O*count1;

        for(int count2 = 0; count2<O+1; count2++){
            Hz=2/(double)O*count2;

            for(int count3 = 0; count3<1; count3++){
                // knowing when to start sampling
                int sampeling_started = 0;
                int samples_taken = 0;
                int max_samples = 30;
                Energy = get_energy_tot(spin);

                char f2[128];
                sprintf(f2, "Data/skyrmions__D_%lf__Hz_%lf__I_%d.json",D, Hz, count3);
                FILE *fp_sk = fopen(f2, "w");

                double betas[] = {0.5, 1.0, 2.0, 4.0};
                int big = sizeof(betas) / sizeof(betas[0]);

                int tot =0;
                for(int s = 0; s<big; s++){
                    beta = betas[s];
                    // for understanding when the energy is stable
                    double E_sum = 0.0;
                    double E_avg_prev = 0.0;
                    int stable_count = 0;
                    int is_stable = 0;


                    for(int count=1; count<M+1; count++){
                        
                        if(count % 1000 == 0 && sampeling_started == 1){
                            Q = get_Q();

                            gaussian_filter_CD(0.8);

                            desission_peak(tot, f2);

                            samples_taken ++;

                            if (samples_taken >= max_samples){
                                break;
                            }
                        }
                        n1 = floor(N*dsfmt_genrand());
                        n2 = floor(N*dsfmt_genrand());
                        accepted += change_particle();
                        tot+=1;

                        E_sum += Energy;

                        if (count % window_size == 0){

                            double E_avg = E_sum / window_size;

                            tolerance = 0.15 * fabs(E_avg);

                            if (count > window_size){ // start at the second windwo 

                                if (fabs(E_avg - E_avg_prev)< tolerance){

                                    stable_count ++;

                                }

                                else stable_count = 0;

                            }

                            E_avg_prev = E_avg; // reset

                            E_sum = 0;

                        }

                        if (stable_count >= number_windows && is_stable == 0){
                            // printf("Stable reached at beta=%f count=%d\n", beta, count);

                            is_stable = 1;

                            if (s < big - 1){
                                // printf("Moving to next beta\n");
                                break;
                            }
                            else{
                                // printf("Starting sampling\n");
                                sampeling_started = 1;
                            }
                        }

                        
                    }
                    
                }

                fclose(fp_sk);
            }
        }
    }

    WriteState2File();

    return 0;
}

/* annealing algorimte to solve a sudoku. we start with a puzzle from wikipedia and fill the rest at random.
We calcualte the amount of conflicts and use that to calcualte the 'energy' and perform a monte carlo algoritm*/


#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"
#include <stdbool.h>


#define N 9  // size of the grid
int M_steps = 4000000; // amount of simulation steps
int M_check = 10000; // check interval
int grid[N][N]; //the grid that is going to be changed
bool fixed[N][N]; // the positions of where the fixed  
double beta = 0.5; //starting temp
double beta_change = 1.01; // increase 5 procent

int puzzle[N][N] = { // puzzle form wikipedia: https://en.wikipedia.org/wiki/Sudoku
    {5,3,0,0,7,0,0,0,0},
    {6,0,0,1,9,5,0,0,0},
    {0,9,8,0,0,0,0,6,0},
    {8,0,0,0,6,0,0,0,3},
    {4,0,0,8,0,3,0,0,1},
    {7,0,0,0,2,0,0,0,6},
    {0,6,0,0,7,0,2,8,0},
    {0,0,0,4,1,9,0,0,5},
    {0,0,0,0,8,0,0,7,9},
};

void fill_grid(void){
   
    for (int r = 0; r<N; r++){
        for(int c = 0; c<N; c++){
            if (puzzle[r][c] != 0){
                grid[r][c]= puzzle[r][c];
                fixed[r][c]= true; // indicat the sqare is fixed
            }
            else{
                grid[r][c]= (int)(9 * dsfmt_genrand())+1; // randomly fill the open square
                fixed[r][c]= false; // indicate the square is open
            }
        }
    }
}

void print_grid(void) {
    printf("Current Sudoku grid:\n");
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            printf("%d ", grid[r][c]);
        }
        printf("\n");
    }
}

void print_fixed(void) {
    printf("fixed squares:\n");
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            printf("%d ", fixed[r][c]);
        }
        printf("\n");
    }
}

int compute_energy(int grid[N][N]){
    /* for the energy we use the amount of conflicts (or mistakes) we made in our sudoku 
    we count how often each number is in a row, collom and square
    afterards we count the amount of conflicts (where count is larger than 1 )*/
   
    int count[10] = {0};
    int conflicts = 0;

    /*checking rows*/
    for (int r = 0; r<N; r++){ // loop over rows
        for( int n = 1; n<= N; n++)count[n]= 0; // set count to zero 
        for(int c = 0; c< N ; c++){  // count over the colom
                count[grid[r][c]]++;  
            }
          for (int n = 1; n<=N; n++ ){ // add the conflict
            if (count[n] >1) conflicts += count[n]-1;
          }
    }

    /*checking coloms*/
    for (int c = 0; c<N; c++){
        for( int n = 1; n<= N; n++) count[n]= 0;
        for(int r = 0; r< N ; r++){
                count[grid[r][c]]++; 
        }
        for (int n = 1; n<=N; n++ ){
            if (count[n] >1) conflicts += count[n]-1;
        }
    }
  

    /*checking blocks*/
    for (int start_r = 0; start_r< N; start_r+= 3){ // loop over the corners of the blocks
        for (int start_c = 0; start_c< N; start_c+= 3){
                
            for( int n = 1; n<= N; n++) count[n]= 0; // set count to 0 for each block
                
            for(int r_b = 0; r_b< 3; r_b++){ // loop over squares int the block
                for(int c_b = 0; c_b< 3; c_b++){

                    int r = start_r+ r_b; // add block ofset
                    int c = start_c+ c_b; // add block ofset

                    count[grid[r][c]]++; // add counter inside the block
                    }
                }
            for (int n = 1; n <= N; n++) {
            if (count[n] > 1) conflicts += count[n] - 1; // add conflicts of each block
            }  
            }
        }

    return conflicts;
}

int change_value (){

    int trail_r = (int)(9 * dsfmt_genrand()); // pick a random square 
    int trail_c = (int)(9 * dsfmt_genrand());

    if (fixed[trail_r][trail_c]) return 0; // reject the change if it is a fixed value of the puzzle

    int old_energy = compute_energy(grid);

    int trail_grid[N][N]; // make a trail grid
    for(int r = 0; r< N; r++){
        for (int c = 0; c< N; c++){
            trail_grid[r][c] = grid[r][c];
        }
    }
    int trail_number= (int)(9 * dsfmt_genrand())+1;

    trail_grid[trail_r][trail_c] = trail_number; // randomly change the one particle

    int new_energy = compute_energy(trail_grid);
    double acc = exp(- beta * (new_energy- old_energy)); // acceptence rule 
   
    if (compute_energy(trail_grid) == 0){ // always accept the right solution
        grid[trail_r][trail_c]= trail_number;
        return 1; 
    }
    else if (dsfmt_genrand()< acc){ // accept trail move
        grid[trail_r][trail_c]= trail_number;
        return 1; 
    }

    return 0; 
}





int main(){

    dsfmt_seed(time(NULL));
    fill_grid();
    print_grid();
    print_fixed();
    printf("amount of conficlts:\n");
    printf("%d \n", compute_energy(grid) );

    float accepted = 0;
    int solved = 0;

    for (int step = 0; step< M_steps; step++){
       
        for(int i = 0; i<N*N; i++){
           accepted += change_value();

           if (compute_energy(grid) == 0){ // check if it is solved
            printf("Solved at step %d\n", step);
            solved = 1;
            break;
            }
        }

        /* terminate if solved */
        if (solved == 1) break;


        if (step % M_check == 0){
            printf("step %d amount of conficlts: %d, beta %.3f acceptence rate: %.3f \n", step, compute_energy(grid), beta, accepted/(M_check* N * N) );
            beta *= beta_change;
            accepted = 0;
        }
    }
    print_grid();
    printf("amount of conficlts:\n");
    printf("%d \n", compute_energy(grid) );
    return 1;
}

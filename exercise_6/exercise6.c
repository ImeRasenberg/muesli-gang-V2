

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define NDIM 3
#define N 1000

/* Initialization variables */
const double diameter = 1.0;

const char* init_filename = "final.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];

/*variables for g(r)*/
#define NBINS 300 // amount of bins
double dr_bin; // size of bins
int nhis[NBINS]; 
double g[NBINS];


/* Functions */

void calcualte_gr(void){

    dr_bin = (0.5* box[0])/NBINS; // Under nearest image max distence is half of the box.
    for (int b = 0; b<NBINS; b++) nhis[b]= 0 ; // Set the histto 0. 
       
    /*fill the hist */
    for (int i = 1; i < n_particles; i++){
        for(int j = 0; j<i; j++){
            double dist2 = 0.0;
            for(int d = 0; d<NDIM; d++){
                double dx = (r[i][d]-r[j][d]);
                dx = dx - box[d] * round(dx /box[d]);
                dist2 += dx *dx;
            }
        double dist = sqrt(dist2);
        int bin = (int)(dist/dr_bin); // Find the bin of the distance.
        nhis[bin]++; // Add one to the bin.


        }
    }
    double rho = n_particles/(box[0]*box[1]* box[2] );
    for(int d = 0; d<NBINS; d++){
       double dis = d*dr_bin; 
       double nid = (4*M_PI*rho/3)*(pow(dis + dr_bin, 3.0) -pow(dis, 3.0) ); // Calcualte ideal gass.
       g[d] =2* nhis[d]/((double)(n_particles)*nid);

    }
}


void read_data(void){
   FILE* coords;
   coords = fopen(init_filename, "r");

   fscanf(coords, "%i\n", &n_particles); // read in number of particles

   for (int i = 0; i < NDIM; i++) {// read in dimension of box 
    double L, R;
    fscanf(coords, "%lf %lf", &L, &R); 
    box[i] = R-L;
    }

    for (int j = 0; j < n_particles; j++) { // read in partilces pos
    fscanf(coords, "%lf %lf %lf %*lf ", &r[j][0], &r[j][1], &r[j][2]);
    }
    fclose(coords);
    /* print out some of the read data as a sanity check*/
    printf("number of particls %d\n", n_particles);
    printf("box size %lf\n", box[0]);
    for (int i = 0; i<3; i++){
        printf("partilce %d\t ,positions: %lf\t%lf\t%lf\n",i, r[i][0], r[i][1], r[i][2] );
    }

}

void write_gr(void){
    FILE* fp = fopen("gr.dat", "w" );
    fprintf(fp, "r   g(r)\n");
    for (int b = 0; b<NBINS; b++){
        double r = (b+ 0.5) * dr_bin;
        fprintf(fp, "%lf\t%lf\n", r, g[b]);
    }
    fclose(fp);
    printf("g(r) written to gr.dat");
    
}
void write_data(){
    char buffer[128];
    sprintf(buffer, "final.dat");
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}


int main(int argc, char* argv[]){
    
    assert(diameter > 0.0);
    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    calcualte_gr(); // calcualte g(r)
    write_gr(); // write it into a file 

    return 0;
}

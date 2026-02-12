#include <stdio.h>
#include <math.h>
// in this file we will make the cubic latice

int main(){
    int N = 4; // The number of particles in each dirrection
    float d = 1.0; // the distance between two spheres
    float a = 1.0; // the radius of an sphere

    // Make a file where we can save the position data
    FILE *print_coords; // inititialises a file variable
    print_coords = fopen("cubic_xyz.dat","w"); // defining the file variable to be the opening of some file cubic.xyz

    // Let us print some initial coordinates
    fprintf(print_coords, "%i\n", N*N*N); // the total number of particles
    fprintf(print_coords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the x direction
    fprintf(print_coords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the y direction
    fprintf(print_coords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the z direction

    // we first initialise the particle possision saving arrays
    float x[N*N*N], y[N*N*N], z[N*N*N], r[N*N*N];


    // now we start generating particle possitions and radiuses
    int n = 0; // this is our counting variable, it wil index which particle we will consider

    // sweeping over the N_x particles
    for(int i=0; i<N; i++){
        // sweeping over the N_y particles
        for(int j=0; j<N; j++){
            // sweeping over the N_z particles
            for(int k=0; k<N; k++){
                // generating the possition for i,j,k latice cite, also the radius of the particle
                x[n]= (i+0.5)*d;
                y[n]= (j)+0.5*d; 
                z[n]= (k+0.5)*d;
                r[n]= a;

                // saving the x,y,z possition and radius of the particle
                fprintf(print_coords, "%lf\t%lf\t%lf\t%lf\n", x[n], y[n],z[n],r[n]);

                n++;
            }
        }
    }


    fclose(print_coords);
    return 0;
}
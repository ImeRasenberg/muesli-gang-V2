#include <stdio.h>
#include <math.h>

int main(){
    int NDIM = 3; //the number of dimmentions reading for reading out the files (you cannot change this to switch to 2D because you particles will overlap????)
    int N; //the number of particles that are read will be read from the files

    // monte carlo simmulation parameters
    int mc_steps = 3;
    
    // the file that will be considerd
    char init_filename[] = "cubic_xyz.dat";

    // degining the file
     FILE *read_cords;
     read_cords = fopen(init_filename, "r");

    // reading the first line to get the number of of particles that exist in the file (why is the exersise so weird???)
    fscanf(read_cords, "%i\n", &N); // the total number of particles
    printf("%i\n", N);

    //defining the space where the read particles will be gotten form the files
    float r[N][NDIM]; // all the position vectors of all the particles
    float box[2][NDIM]; // the size of the box 

    

    // fscanf(read_cords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the x direction
    // fscanf(read_cords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the y direction
    // fscanf(read_cords, "%lf\t%lf\n", -0.0, 1.0*d*N); // The ocupied space in the z direction


    // fscanf(read_cords,),


    fclose(read_cords);
    return 0;
}
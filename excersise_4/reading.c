#include <stdio.h>
#include <math.h>
int read_data(){
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
    float box[2][NDIM]; // the size of the box (to make the particles fit inside the box poroperly 2 points are defined for me)

    // apearantly you cannot read long floats you just have to read them as floats???
    // fscanf(read_cords, "%f\t%f\n", &box[0][0], &box[1][0]);
    // printf("%lf\t%lf\n", box[0][0], box[1][0]);
    // fscanf(read_cords, "%f\t%f\n", &box[0][1], &box[1][1]);
    // printf("%lf\t%lf\n", box[0][1], box[1][1]);
    // fscanf(read_cords, "%f\t%f\n", &box[0][2], &box[1][2]);
    // printf("%lf\t%lf\n", box[0][2], box[1][2]);

    // lets turn the above into a loop because i want to
    for(int i = 0; i<NDIM; i++){
        // This reads the Min into box[0][i] and Max into box[1][i]
        fscanf(read_cords, "%f %f", &box[0][i], &box[1][i]);
        printf("%f %f\n", box[0][i], box[1][i]);
    }

    // now that we have arived at the paricles lets be happy
    for(int i = 0; i<N; i++){
        fscanf(read_cords, "%f %f %f", &r[i][0], &r[i][1], &r[i][2]);
    }


    fclose(read_cords);
    return 0;

}


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
    float box[2][NDIM]; // the size of the box (to make the particles fit inside the box poroperly 2 points are defined for me)

    // apearantly you cannot read long floats you just have to read them as floats???
    // fscanf(read_cords, "%f\t%f\n", &box[0][0], &box[1][0]);
    // printf("%lf\t%lf\n", box[0][0], box[1][0]);
    // fscanf(read_cords, "%f\t%f\n", &box[0][1], &box[1][1]);
    // printf("%lf\t%lf\n", box[0][1], box[1][1]);
    // fscanf(read_cords, "%f\t%f\n", &box[0][2], &box[1][2]);
    // printf("%lf\t%lf\n", box[0][2], box[1][2]);

    // lets turn the above into a loop because i want to
    for(int i = 0; i<NDIM; i++){
        // This reads the Min into box[0][i] and Max into box[1][i]
        fscanf(read_cords, "%f %f", &box[0][i], &box[1][i]);
        printf("%f %f\n", box[0][i], box[1][i]);
    }

    // now that we have arived at the paricles lets be happy
    for(int i = 0; i<N; i++){
        fscanf(read_cords, "%f %f %f", &r[i][0], &r[i][1], &r[i][2]);
    }


    fclose(read_cords);
    return 0;
}
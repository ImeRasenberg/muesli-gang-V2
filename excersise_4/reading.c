#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
// // //example found online

#include <stdio.h>

typedef struct {
    float **A; // makes matrix A_{i,j} each star gives an index
    float b;
} Result;

Result create_data() {

    Result r;

    r.A[0][0] = 1.0;
    r.A[0][1] = 2.0;
    r.A[1][0] = 3.0;
    r.A[1][1] = 4.0;

    r.b = 5.0;

    return r;
}

int main() {

    Result data = create_data();

    printf("Matrix:\n");
    printf("%f %f\n", data.A[0][0], data.A[0][1]);
    printf("%f %f\n", data.A[1][0], data.A[1][1]);

    printf("scalar:\n")
    printf("b = %f\n", data.b);

    return 0;
}

*/

/*
here we see we can make a new file structure and define a function to output that structure

apearently its better to define rows and colombs in one pointer....? why?

then we define a function that puts data into the struture

finally we run it in the main function

this is what we will do for this code
*/
// #defining the datastructure of interest
typedef struct
{
int N; // getting the number of particles from the file
float (*box)[2]; //the size of the box, one direction doesnt yet have a defined size
float (*r)[3]; // same idea here but with the possition of the particle
float *size; // the number of particles
} Loaded_Data;

Loaded_Data load_data( char *init_filename){
    Loaded_Data Loaded_Data;
    int NDIM=3;
    // degining the file
     FILE *read_cords;
     read_cords = fopen(init_filename, "r");

    // reading the first line to get the number of of particles that exist in the file (why is the exersise so weird???)
    fscanf(read_cords, "%i\n", &Loaded_Data.N); // the total number of particles
    printf("%i\n", Loaded_Data.N);

    //defining the space where the read particles will be gotten form the files
    // float r[Loaded_Data.N][NDIM]; 
    // float box[2][NDIM];     // float size[Loaded_Data.N]; 

    // making sure that the size will be correctly degined instead of having to asign it before hand
    // malloc is the memmory allocation commman which is wat we need to have exact size matrixes, only this satisfies me
    Loaded_Data.box = malloc(NDIM * sizeof * Loaded_Data.box); // the size of the box (to make the particles fit inside the box poroperly 2 points are defined for me)
    Loaded_Data.r = malloc(Loaded_Data.N * sizeof * Loaded_Data.r); // all the position vectors of all the particles
    Loaded_Data.size = malloc(Loaded_Data.N * sizeof * Loaded_Data.size); //The size of all particles

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
        fscanf(read_cords, "%f %f", &Loaded_Data.box[0][i], &Loaded_Data.box[1][i]);
        printf("%f %f\n", Loaded_Data.box[0][i], Loaded_Data.box[1][i]);
    }

    // now that we have arived at the paricles lets be happy
    for(int i = 0; i<Loaded_Data.N; i++){
        fscanf(read_cords, "%f %f %f %f", &Loaded_Data.r[i][0], &Loaded_Data.r[i][1], &Loaded_Data.r[i][2], &Loaded_Data.size[i]);
    }

    fclose(read_cords);
    return Loaded_Data;
}


typedef struct
{
int Index;
float displacement
} displacement;

displacement move_particle(float Delta){
    displacement d;


    
    return d;
    
}

int main(){
    int NDIM = 3; //the number of dimmentions reading for reading out the files (you cannot change this to switch to 2D because you particles will overlap????)

    // the file that will be considerd
    char *init_filename= "cubic_xyz.dat";

    Loaded_Data Loaded_Data = load_data(init_filename);



    // monte carlo simmulation parameters
    int mc_steps = 3;
    
    return 0;
}
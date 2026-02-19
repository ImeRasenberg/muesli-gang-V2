#include <stdio.h>
#include <stdlib.h>

 
#include <math.h>

#include <time.h> // we use time as a random seed
#include "../downloads/mt19937.h"

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
    // printf("%i\n", Loaded_Data.N);

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
        fscanf(read_cords, "%f %f", &Loaded_Data.box[i][0], &Loaded_Data.box[i][1]);
        // printf("%f %f\n", Loaded_Data.box[0][i], Loaded_Data.box[1][i]);
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
int index;
float d_r[3];
float l;
int disp;
} displacement;

int check_particle_overlap(Loaded_Data l, displacement d) {
    // getting the selected particle and adding the gotten displacement to it
    float *p = l.r[d.index];
    // p[0] += d.d_x;
    // p[1] += d.d_y;
    // p[2] += d.d_z;

    // checking over all particles
    for (int i = 0; i < l.N; i++) {
        // if we consider the same particle the distance is always smaller then the combination of the two raduses
        if (i == d.index) {
            continue;
        }
        // getting the postion of the particle
        float *p_c = l.r[i];

        // setting the distance between this particle and the changed particle to 0
        float dist_sq = 0.0;
        // updating this dastance to be acurate
        for (int j = 0; j < 3; j++) {
            float diff = p[j] + d.d_r[j] - p_c[j];
            float length = abs(l.box[j][1] - l.box[j][0]);
            if (diff>0.5*length){
                diff -= length;
                
            }
            else if (diff<-0.5*length){
                diff += length;
            }

            dist_sq += diff * diff;
        }

        // making degining the minimum distance between two
        float sum_raduses = 0.5 * (l.size[i] + l.size[d.index]);
        

        if (dist_sq < sum_raduses*sum_raduses) {
            // printf("%f\t%f\t%f\nf",d.d_r[0],d.d_r[1],d.d_r[2]);
            return 1; // Overlap detected
        }
    }
    // printf("there is no overlap\n");
    return 0; // No overlaps found
}



displacement move_particle(float Delta,Loaded_Data Loaded_Data){
    displacement d;

    d.index =floor(dsfmt_genrand()*Loaded_Data.N);
    d.d_r[0] = (dsfmt_genrand()-0.5);
    d.d_r[1] = (dsfmt_genrand()-0.5);
    d.d_r[2] = (dsfmt_genrand()-0.5);

    d.l = (dsfmt_genrand()-0.5)*Delta+0.0001; //no devision by 0

    float length = sqrt(d.d_r[0]*d.d_r[0]+d.d_r[1]*d.d_r[1]+d.d_r[2]*d.d_r[2]);

    d.d_r[0] = d.d_r[0]/length*d.l;
    d.d_r[1] = d.d_r[1]/length*d.l;
    d.d_r[2] = d.d_r[2]/length*d.l;


    d.disp =  check_particle_overlap(Loaded_Data, d);


    // printf("%i\n",displacement.disp);
    if(d.disp == 1){
        // printf("overlap found no displacement\n");
        return d; 
    }
    else if (d.disp == 0)
    {
        // printf("from: \t%f\t%f\t%f\n",Loaded_Data.r[displacement.index][0],
        //     Loaded_Data.r[displacement.index][1],Loaded_Data.r[displacement.index][2]);
        Loaded_Data.r[d.index][0] += d.d_r[0];
        Loaded_Data.r[d.index][1] += d.d_r[1];
        Loaded_Data.r[d.index][2] += d.d_r[2];

        for(int l=0; l<3;l++){
            float box_len = abs(Loaded_Data.box[l][0]-Loaded_Data.box[l][1]) ;

            if(Loaded_Data.r[d.index][l]<Loaded_Data.box[l][0]){
                Loaded_Data.r[d.index][l]+= box_len;
            }
            if(Loaded_Data.r[d.index][l]>Loaded_Data.box[l][1]){
                Loaded_Data.r[d.index][l] -= box_len;
            }
        }
        // printf("to: \t%f\t%f\t%f\n",Loaded_Data.r[displacement.index][0],
        // Loaded_Data.r[displacement.index][1],Loaded_Data.r[displacement.index][2]);
        // printf("done a displacement\n");

        
    }


    return d;
}



void write_to_file(Loaded_Data l){

    FILE *print_coords; // inititialises a file variable
    char *new_name = "NVT_output.dat";
    print_coords = fopen(new_name,"w");

    fprintf(print_coords, "%i\n", l.N); // the total number of particles
    fprintf(print_coords, "%lf\t%lf\n", l.box[0][0], l.box[0][1]); // The ocupied space in the x direction
    fprintf(print_coords, "%lf\t%lf\n", l.box[1][0], l.box[1][1]); // The ocupied space in the y direction
    fprintf(print_coords, "%lf\t%lf\n", l.box[2][0], l.box[2][1]); // The ocupied space in the z direction

    for(int i=0;i<l.N;i++){

        fprintf(print_coords, "%lf\t%lf\t%lf\t%lf\n", l.r[i][0], l.r[i][1],l.r[i][2],l.size[i]);
    }

    fclose(print_coords);
}




int main(){
    dsfmt_seed(time(NULL)); //setting the seed for the random displacement

    int NDIM = 3; //the number of dimmentions reading for reading out the files (you cannot change this to switch to 2D because you particles will overlap????)

    int succes_count=0 ;
    int mc_steps = 100000;

    int dV_m = 0.1;
    // the file that will be considerd
    char *init_filename= "FCC_xyz.dat";

    Loaded_Data Loaded_Data = load_data(init_filename);
    printf("starting displacement loop\n");
    for (int k=0; k< mc_steps; k++){

        // printf(" run %i\n",k);
        displacement displacement = move_particle(0.1, Loaded_Data);
        
        if (displacement.disp == 0){
            succes_count+=1;
        }
        
    }
    printf("finnished displacement loop\n");

    printf("fracction succes: %i/%i\n",succes_count,mc_steps);
    write_to_file(Loaded_Data);
    return 0;
}
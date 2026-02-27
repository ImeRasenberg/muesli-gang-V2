#include <stdio.h>
#include <math.h>

int main(){
    
    int N=6;
    float d=1.0;
    float a=sqrt(2.0)*d;

    FILE *print_coords;
    print_coords=fopen("fcc.xyz","w");
    fprintf(print_coords, "%i\n", 4*N*N*N);
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a);
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a);
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a);

    float x1[N*N*N], y1[N*N*N], z1[N*N*N];
    float x2[N*N*N], y2[N*N*N], z2[N*N*N];
    float x3[N*N*N], y3[N*N*N], z3[N*N*N];
    float x4[N*N*N], y4[N*N*N], z4[N*N*N];

    int n=0;

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){

                x1[n]=i*a;
                y1[n]=j*a;
                z1[n]=k*a;
                fprintf(print_coords,
                "%lf\t%lf\t%lf\t%lf\n",
                x1[n],y1[n],z1[n],d);

                x2[n]=i*a+a/2;
                y2[n]=j*a+a/2;
                z2[n]=k*a;
                fprintf(print_coords,
                "%lf\t%lf\t%lf\t%lf\n",
                x2[n],y2[n],z2[n],d);

                x3[n]=i*a+a/2;
                y3[n]=j*a;
                z3[n]=k*a+a/2;
                fprintf(print_coords,
                "%lf\t%lf\t%lf\t%lf\n",
                x3[n],y3[n],z3[n],d);
        
                x4[n]=i*a;
                y4[n]=j*a+a/2;
                z4[n]=k*a+a/2;
                fprintf(print_coords,
                "%lf\t%lf\t%lf\t%lf\n",
                x4[n],y4[n],z4[n],d);

                n++;

            }
        }
    }
    fclose(print_coords);
}

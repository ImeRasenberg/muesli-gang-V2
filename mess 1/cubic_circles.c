#include <stdio.h>
#include <math.h>
int main(){
    int N=5;
    float a=1.0;
    float d=1.0;

    FILE *print_coords;
    print_coords = fopen("cubic.xyz","w");

    fprintf(print_coords, "%i\n", N*N*N);
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );

    float x[N*N*N];
    float y[N*N*N];
    float z[N*N*N];

    int n=0;
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
            
                x[n]=i*a;
                y[n]=j*a;
                z[n]=k*a;

                fprintf(print_coords, "%lf\t%lf\t%lf\t%lf\t\n", x[n],y[n],z[n],d);

                n++;
            }
        }
    }
    fclose(print_coords);

    return 0;
}
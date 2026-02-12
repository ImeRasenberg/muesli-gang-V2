#include <stdio.h>
#include <math.h>

int main(){

    int N=4;
    float d=1.0;
    float a = sqrt(2.0)*d;

    FILE *print_coords;
    print_coords = fopen("fcc.xyz","w");
    fprintf(print_coords, "%i\n", 4*N*N*N);
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );
    fprintf(print_coords, "%lf\t%lf\n", 0.0, (N)*a );

    float x1[N,N,N], y1[N,N,N], y1[N,N,N];
}
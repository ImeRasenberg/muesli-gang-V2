#include <stdio.h>
#include <math.h> 
#include <time.h> // we use time as a random seed
// #include "downloads/mt19937.h"
#include "../downloads/mt19937.h"


int main(){
    dsfmt_seed(time(NULL));

    int inside = 0;
    int N = 1E6;

    for( int i=0; i < N; i++){
        double x = dsfmt_genrand();
        double y = dsfmt_genrand();
        double r2 = x*x + y*y;

        if(r2<=1.0*1.0) inside++;
    }

    double pi = 4.0 *((double)(inside))/((double)(N))/1.0;

    printf("%lf", pi);

    return 0;
}


//https://webspace.science.uu.nl/~herme107/viscol/
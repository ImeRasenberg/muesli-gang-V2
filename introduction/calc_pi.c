#include <math.h>
#include <stdio.h>


int main(){
    double x,y,z;
    double pi;

    pi = atan(1)*4;
    // x = exp(atan(1)*4);
    x = exp(pi);

    y = sinh(x);
    z = pow(x,y);


    printf("pi is %f! (only ! isnt the factorial)\n", pi);
    return 0;
}

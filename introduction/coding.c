#include <stdio.h>

int main2(){
    int i;
    for (i=1;i<20;i++)
        if (i%2 == 0) printf("%d is even!\n" , i);
        else printf("%d is odd!\n" , i);
    return 0;
}


int main(){
    int i;

    for  (i=1; i<20; i++)
        {
            if (i % 2 ==0)
            {
                printf("%d is even!\n" , i);
            }
            else
            {
                printf("%d is odd!\n" , i);
            }
        }
        main2();
    
    return 0;
}

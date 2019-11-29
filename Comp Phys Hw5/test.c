#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <stdbool.h> //standard boolean library
#include <string.h>  //string routines
#include <math.h>    //math routines

void main()
{
    int i;
    float f[101];
    float fs[101];
    float ft[101];
    float del[101];
    float dif = 0;

    f[0] = 0;
    f[100] = 0;
    fs[0] = f[0];
    fs[100] = f[100];

    for(i = 1; i <= 99; i++)
    {
        f[i] = 0;
    }

    do
    {
        dif = 0;

        for(i = 1; i <= 99; i++)
        {
            fs[i] = (f[i - 1] + f[i + 1]) / 2 + 4.9 * 0.01;
            del[i] = fs[i] - f[i];
        }

        for(i = 1; i <= 99; i++)
        {
            if(dif < fabs(del[i]))
            {
                dif = fabs(del[i]);
            }
        }
        //printf("%f", dif);

        for(i = 0; i <= 100; i++)
        {
            ft[i] = fs[i];
            fs[i] = f[i];
            f[i] = ft[i];
        }
        //printf("%9.8f\n", f[50]);
    } while (dif > 0.000001);

    for(i = 0; i <= 100; i++)
    {
        printf("%f\n", f[i]);
    }
    //printf("%f", f[50]);
}

/*
root_finding.h

Header file for root finding methods

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <string.h>  //string routines
#include <math.h>    //math routines
*/

char equation_name[64];

/*
Binary search method

To calculate the root of equation f(x) = 0, start with two points x_1 and x_2, satisfying f(x_1) * f(x_2) <= 0
Then bisect the range and find if the midpoint x' satisfies f(x_1) * f(x') <= 0 or f(x') * f(x_2) <= 0
Keep x_1 or x_2 that satisfies the relation above and renew that half interval.
Repeat these steps until reaching the target accuracy |x_2 - x_1| < accuracy

This function returns the process of root finding and the final solution.
*/
float bin_search(float x_1, float x_2, float (*func_ptr)(float x), float accuracy)
{
    FILE *fp;
    char filename[64];
    float x;            //midpoint of the interval

    sprintf(filename, "bin_search_of_%s=0.txt", equation_name);
    fp = fopen(filename, "w");
    fprintf(fp, "Solution to %s=0 using binary search\n\n", equation_name);
    fprintf(fp, "x_1 = %7.6e, x_2 = %7.6e, accuracy = %7.6e\n\n", x_1, x_2, accuracy);
    fprintf(fp, "      x        %s      \n", equation_name);

    printf("%7.6e %7.6e\n", x_1, (*func_ptr)(x_1));
    printf("%7.6e %7.6e\n", x_2, (*func_ptr)(x_2));
    fprintf(fp, "%7.6e %7.6e\n", x_1, (*func_ptr)(x_1));
    fprintf(fp, "%7.6e %7.6e\n", x_2, (*func_ptr)(x_2));

    do
    {
        x = (x_1 + x_2) / 2;
        printf("%7.6e %7.6e\n", x, (*func_ptr)(x));
        fprintf(fp, "%7.6e %7.6e\n", x, (*func_ptr)(x));

        if((*func_ptr)(x_1) * (*func_ptr)(x) <= 0)  //if f(x_1) * f(x') <= 0, it means there must be a root within [x_1, x], or the root will be within [x, x_2]
            x_2 = x;
        else
            x_1 = x;
    } while(fabs(x_2 - x_1) >= accuracy);

    fprintf(fp, "\nx = %7.6e", (x_1 + x_2) / 2);

    return((x_1 + x_2) / 2);

    fclose(fp);
}

/*
relaxation.h

Solution of nonlinear equations using relaxation methods.

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <string.h>  //string routines
#include <math.h>    //math routines
*/

char equation_name[64];

/*
Solving one dimensional nonlinear equation x = f(x) using relaxation method : new x' = f(x) until convergence.
Error \epsilon \approx (x - x') / (1 - 1 / f'(x))
It is required |f'(x*)| < 1 to converge.
*/
float relaxation(float x_0, float (*func_ptr)(float x), float accuracy)
{
    float x;            //old estimate of x
    float x_prime;      //new estimate of x
    float epsilon;      //approximate error
    int iter;           //number of iterations
    FILE *fp;
    char filename[64];

    sprintf(filename, "relax_of_%s.txt", equation_name);
    fp = fopen(filename, "w");
    fprintf(fp, "Solution to %s using relaxation method\n\n");
    fprintf(fp, "x_0 = %7.6e\n\n", x_0);
    fprintf(fp, "iter     x        \n");

    iter = 1;
    x = (*func_ptr)(x_0);
    epsilon = (x - x_0) / (1 - (x - x_0) / ((*func_ptr)(x) - x));  //x = f(x_0)
    printf("  i         x          \n");
    printf("%3d %7.6e\n", iter, x);
    fprintf(fp, "%3d %7.6e\n", iter, x);
    iter++;

    for(; fabs(epsilon) >= fabs(accuracy); iter++)
    {
        x_prime = (*func_ptr)(x);
        epsilon = (x_prime - x) / (1 - (x_prime - x) / ((*func_ptr)(x_prime) - x_prime));  //x' = f(x)
        x = x_prime;
        printf("%3d %7.6e\n", iter, x);
        fprintf(fp, "%3d %7.6e\n", iter, x);
    }

    fprintf(fp, "\nx = %7.6e\n\n", x);

    fclose(fp);

    return x;
}

/*
Solving one dimensional nonlinear equation x = f(x) using overrelaxation method

x' = (1 + \omega)f(x) - \omega * x

\epsilon \approx (x - x') / (1 - 1 / [(1 + \omega) * f'(x) - \omega])

This function will give the result for some different \omega's and also return the number of iterations
for every \omega to reach the target accuracy
*/
float overrelaxation(float x_0, float (*func_ptr)(float x), float accuracy)
{
    float x;            //old estimate of x
    float x_prime;      //new estimate of x
    float epsilon;      //approximate error
    int iter;           //number of iterations
    float omega;        //overrelaxation parameter
    FILE *fp;
    char filename[64];

    sprintf(filename, "overrelax_of_%s.txt", equation_name);
    fp = fopen(filename, "w");
    fprintf(fp, "Solution to %s using overrelaxation method\n\n");
    fprintf(fp, "x_0 = %7.6e\n\n", x_0);
    fprintf(fp, "iter     \\omega        x        \n");
    printf("iter     \\omega        x        \n");

    for(omega = 0; omega < 1; omega += 0.01)
    {
        iter = 1;
        x = (*func_ptr)(x_0) * (1 + omega) - omega * x_0;
        epsilon = (x - x_0) / (1 - 1 / ((1 + omega) * ((*func_ptr)(x) - x) / (x - x_0) - omega));  //x = f(x_0)
        iter++;

        for(; fabs(epsilon) >= fabs(accuracy); iter++)
        {
            x_prime = (*func_ptr)(x) * (1 + omega) - omega * x;
            epsilon = (x_prime - x) / (1 - 1 / ((1 + omega) * ((*func_ptr)(x_prime) - x_prime) / (x_prime - x) - omega));  //x' = f(x)
            x = x_prime;
        }

        printf("%3d %7.6e %7.6e\n", iter - 1,omega, x);
        fprintf(fp, "%3d %7.6e %7.6e\n", iter - 1,omega, x);
    }

    fclose(fp);

    return x;
}

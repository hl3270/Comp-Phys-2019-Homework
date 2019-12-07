/*
descent.h

Descent method to find local minimums of function, e.g. gradient descent, steepest descent

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routine
*/

char function_name[64];
char variable_name[5][13];
float eps_tol;

/*
Gradient descent method in multi dimensions

To find the minimum of f(*x), start from *x_0

x_n+1[i] = x_n[i] - \gamma * df(*x_n) / dx[i]
*/
float *gradient_descent(int dim, float (*func_ptr)(float *x), float *x_0, float gamma, float eps, int max_loop)
{
    float *x;      //independent variables
    float *x_temp; // x template for derivative calculation
    float *x_new;  //new x vector storage
    float max_dif; //maximum difference of x[i] between two steps
    float h;       //numerical derivative step size h
    float df;      //derivative of f(*x)
    float df_temp; //template of derivative to compare
    FILE *fp;
    char filename[64];

    sprintf(filename, "grad_descent_of_%s.txt", function_name);
    fp = fopen(filename, "w");
    printf("Gradient descent of %s\n\n", function_name);
    fprintf(fp, "Gradient descent of %s\n\n", function_name);
    printf("gamma = %7.6e, epsilon = %7.6e\n\n", gamma, eps);
    fprintf(fp, "gamma = %7.6e, epsilon = %7.6e\n\n", gamma, eps);
    for(int i = 0; i < 1; i++)
    {
        printf("%s ", variable_name[i]);
        fprintf(fp, "%s ", variable_name[i]);
    }
    printf("\n");
    fprintf(fp, "\n");
    for(int i = 0; i <= dim - 1; i++)
    {
        printf("%7.6e ", x_0[i]);
        fprintf(fp, "%7.6e ", x_0[i]);
    }
    printf("%7.6e\n", (*func_ptr)(x_0));
    fprintf(fp, "%7.6e\n", (*func_ptr)(x_0));

    x = (float *) malloc(dim * sizeof(float));
    x_temp = (float *) malloc(dim * sizeof(float));
    x_new = (float *) malloc(dim * sizeof(float));

    for(int i = 0; i < dim; i++)       //initialize
    {
        x[i] = x_0[i];
        x_temp[i] = x_0[i];
        x_new[i] = x_0[i];
    }

    int j = 0;
    do
    {
        max_dif = 0;
        for(int i = 0; i < dim; i++)   //calculate gradient
        {
            h = 0.05;
            x_temp[i] = x[i] + h;      //numerical derivative using midpoint rule
            df_temp = func_ptr(x_temp) / (2 * h);
            x_temp[i] = x[i] - h;
            df_temp -= func_ptr(x_temp) / (2 * h);
            h /= 2;

            do
            {
                df = df_temp;
                x_temp[i] = x[i] + h;
                df_temp = func_ptr(x_temp) / (2 * h);
                x_temp[i] = x[i] - h;
                df_temp -= func_ptr(x_temp) / (2 * h);
                h /= 2;
            } while (fabs(df_temp - df) >= eps_tol);

            x_temp[i] = x[i];
            x_new[i] -= gamma * df_temp;
            //printf("%f\n", df_temp);

            if(max_dif < fabs(gamma * df_temp))
                max_dif = fabs(gamma * df_temp);
        }

        for(int i = 0; i < dim; i++)   //descent
        {
            x[i] = x_new[i];
            x_temp[i] = x_new[i];
        }

        for(int i = 0; i <= dim - 1; i++)
        {
            printf("%7.6e ", x[i]);
            fprintf(fp, "%7.6e ", x[i]);
        }
        printf("%7.6e\n", (*func_ptr)(x));
        fprintf(fp, "%7.6e\n", (*func_ptr)(x));

        j++;
    } while ((max_dif >= eps) && (j < max_loop));

    printf("\nsteps = %d\n",j);
    fprintf(fp, "\nsteps = %d\n",j);

    return x;

    free(x_temp);
    fclose(fp);
}

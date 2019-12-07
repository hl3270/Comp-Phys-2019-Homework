/*Hw_2-3.c*/

/*
Attached is a data file that contains measurements of the galaxy stellar mass functions from the COSMOS galaxy survey.
The columns are: 1. log M_gal [dex]
                 2. n(M_gal)[1/dex/Volume]
                 3. error in n(M_gal)
Here, "dex" means base-10 log of the stellar mass. Volume here is (Mpc/h)^3.

Measurements of this type are usually described by a "Schechter" function, which for this problem has the form

n(M_gal) = \phi^* * (M_gal / M^*) ^ (\alpha + 1) * exp(-M_gal / M^*) * ln(10)

This function (which has the same units as the column 2 in the data file) has three free parameters :
\phi^*, the amplitude,
M^*, the "characteristic mass scale", where the function changes from being a power-law to having an exponential cutoff,
and \alpha, the low-mass slope of the function.

Write a code that will implement the gradient descent method in multiple dimensions,
using numerical derivatives, to find the minimum of a function.
Test your code on a simple function, such as : f(x, y) = (x - 2) ^ 2 + (y - 2) ^ 2
(Remember that we're still using numerical derivatives, even though the derivatives of this function are known.)

After confirming that your code works(include a plot demonstrating this result
- use your discretion to choose what the plot should show), apply this to the data file provided.
The function you are minimizing is the \chi^2 of model.
Verify that your result is robust by demonstrating that you get the same result
when starting from distinct locations in parameter space (within reason).

Attach the following plots :
(a) \chi^2 as a function of step i.
(b) a comparison of your best fit Schechter function to the data, on a log-log plot.
If you want to attach more plots to help you explain your results, fell free to do so.
*/

#include <stdlib.h>        //standard library
#include <stdio.h>         //standard input/output library
#include <string.h>        //string routines
#include <math.h>          //math routines
#include <malloc.h>        //memory allocation routine
#include <stdbool.h>       //standard boolean routine

#include "descent.h"
#include "data_file.h"

/*Structure of the galaxy stellar mass function from the COSMOS galaxy survey*/
struct smf_cosmos
{
    int n;                       //number of data point

    float *log_stellar_mass;     //log10 of stellar mass
    float *n_m_gal;              //density of galaxy with m_gal
    float *error_n;              //error of n_m_gal
};
struct smf_cosmos smf;

void main()
{
    float func_2(float *x);      //the equation f(x) = 0 pending root finding
    float chi_square_n(float *parameter);
    float *gradient_descent(int dim, float (*func_ptr)(float *), float *x_0, float gamma, float eps, int max_loop);  //gradient descent method

    int dim;
    float *x_0;                  //independent variable
    float *x;
    float (*func_ptr)(float *x); //function pointer
    float gamma;
    float eps;                   //target accuracy
    int max_loop;
    FILE *fp;
    float *para;

    /*Test the gradient descent method with a simple function f(x, y) = (x - 2) ^ 2 + (y - 2) ^ 2*/
    dim = 2;
    x_0 = (float *) malloc(dim * sizeof(float));
    x = (float *) malloc(dim * sizeof(float));
    func_ptr = func_2;
    x_0[0] = 0;
    x_0[1] = 1;
    gamma = 1e-3;
    eps = 1e-7;
    max_loop = 6000;
    eps_tol = 5e-5;
    strcpy(function_name, "(x-2)^2+(y-2)^2");
    strcpy(variable_name[0], "      x       ");
    strcpy(variable_name[1], "      y       ");
    strcpy(variable_name[2], "    f(x,y)    ");
    //x = gradient_descent(dim, func_ptr, x_0, gamma, eps, max_loop);
    //printf("\nx = %7.6e, y = %7.6e, f = %7.6e\n\n", x[0], x[1], (*func_ptr)(x));
    free(x);
    free(x_0);

    printf("\n-------------------------------------------------\n\n\n");

    /*Calculate the coefficient of the Schechter function using gradient descent method*/
    dim = 3;
    func_ptr = chi_square_n;
    gamma = 1e-4;
    eps = 2e-7;
    max_loop = 3000;
    eps_tol = 3e-2;
    strcpy(function_name, "Schechter");
    strcpy(variable_name[0], "   log(phi*)   ");
    strcpy(variable_name[1], "     log(M*)   ");
    strcpy(variable_name[2], "      alpha    ");
    strcpy(variable_name[3], "      chi^2    ");
    strcpy(fn.filename, "smf_cosmos.txt");
    fp = fopen(fn.filename, "r");
    smf.n = row_file(fn.filename);
    smf.log_stellar_mass = (float *) malloc(smf.n * sizeof(float));
    smf.n_m_gal = (float *) malloc(smf.n * sizeof(float));
    smf.error_n = (float *) malloc(smf.n * sizeof(float));
    for(int i = 0; i <= smf.n - 1; i++)
    {
        fscanf(fp, "%f %e %e\n", &smf.log_stellar_mass[i], &smf.n_m_gal[i], &smf.error_n[i]);
    }

    para = (float *) malloc(3 * sizeof(float));
    para[0] = - 5;
    para[1] = 9.5;
    para[2] = - 1.5;
    //printf("%f\n", chi_square_n(para));
    para = gradient_descent(dim, func_ptr, para, gamma, eps, max_loop);
    printf("\nchi^2 = %f\n", chi_square_n(para));
    fclose(fp);
    free(para);
    free(smf.error_n);
    free(smf.log_stellar_mass);
    free(smf.n_m_gal);
}

/*Simple test function f(x, y) = (x - 2) ^ 2 + (y - 2) ^ 2. *x = {x, y}*/
float func_2(float *x)
{
    return(pow(x[0] - 2, 2) + pow(x[1] - 2, 2));
}

/*chi^2 of n. *parameter = {log(\phi*), log(M*), \alpha}*/
float chi_square_n(float *parameter)
{
    float n_m;                       //expected value of n_m_gal
    float chi_sq;                    //chi square

    chi_sq = 0;                      //initialize chi square 0
    for(int i = 0; i <= smf.n - 1; i++)
    {
        n_m = pow(10, parameter[0]) * pow(pow(10, smf.log_stellar_mass[i] - parameter[1]), parameter[2] + 1) * exp(- pow(10, smf.log_stellar_mass[i] - parameter[1]))* log(10);
        chi_sq += pow((n_m - smf.n_m_gal[i]), 2) / pow(smf.error_n[i], 2);
    }

    return chi_sq;
}

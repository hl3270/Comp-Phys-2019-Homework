/*Hw_4-1.c*/

/*
Exercise 8.18 (Computational Physics, Mark Newman)
Oscillating chemical reactions

The Belousov–Zhabotinsky reaction is a chemical oscillator, a cocktail of chemicals which, when heated,
undergoes a series of reactions that cause the chemical concentrations in the mixture to oscillate between two extremes.
You can add an indicator dye to the reaction which changes color depending on the concentrations
and watch the mixture switch back and forth between two different colors for as long as you go on heating the mixture.
Physicist Ilya Prigogine formulated a mathematical model of this type of chemical oscillator,
which he called the "Brusselator" after his home town of Brussels. The equations for the Brusselator are

dx / dt = 1 - (b + 1) * x + a * x ^ 2 * y
dy / dt = b * x - a * x ^ 2 * y

Here x and y represent concentrations of chemicals and a and b are positive constants.

Write a program to solve these equations for the case a = 1, b = 3 with initial conditions x = y = 0, to an accuracy of
at least d = 10^(-10) per unit time in both x and y, using the adaptive Bulirsch–Stoer method described in Section 8.5.6.
Calculate a solution from t = 0 to t = 20, initially using a single time interval of size H = 20.
Allow a maximum of n = 8 modified midpoint steps in an interval before you divide in half and try again.

Make a plot of your solutions for x and y as a function of time, both on the same graph, and have your program add dots
to the curves to show where the boundaries of the time intervals lie. You should find that the points are significantly
closer together in parts of the solution where the variables are changing rapidly.

Hint: The simplest way to perform the calculation is to make use of recursion, as described in Exercise 8.17.
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines

#include "ODE.h"

double a = 1;
double b = 3;

void main()
{
    double fx(double x, double y, double t);
    double fy(double x, double y, double t);
    double *mod_mid_2D(double x, double y, double t, double (*fx_ptr)(double x, double y, double t), double (*fy_ptr)(double x, double y, double t), int n, double H0);
    double *Bulirsch_Stoer_2D(double x, double y, double t, double (*fx_ptr)(double x, double y, double t), double (*fy_ptr)(double x, double y, double t), double H0);

    double x;
    double y;
    double *z;               //output of modified midpoint method
    int *lev;                //the unfinished time section points counted from the right
    double t;
    double H0;
    int i;                   //index of unfinished time points
    double h_i;
    FILE *fp;
    char filename[64];

    x = 0;
    y = 0;
    t = t_0;
    //z = mod_mid_2D(x, y, t, fx, fy, 2, 2);
    //printf("%f o\n\n", z[0]);

    /*Solve the Brusselator with modified midpoint method*/
    sprintf(filename, "BZ_reaction_mod_mid.txt");
    fp = fopen(filename, "w");
    fprintf(fp, "The Belousov–Zhabotinsky reaction with modified midpoint method\n");
    fprintf(fp, "dx / dt = 1 - (b + 1) * x + a * x ^ 2 * y\n");
    fprintf(fp, "dy / dt = b * x - a * x ^ 2 * y\n");
    fprintf(fp, "a = 1, b = 3, x_0 = 0, y_0 = 0, accuracy delta = 1e-10, H = 20\n\n");
    fprintf(fp, "       t            x             y\n");
    fprintf(fp, "%7.6e %7.6e %7.6e\n", t, x, y);

    for(i = 1; i <= 2000; i++)
    {
        z = mod_mid_2D(x, y, t, fx, fy, 8, 0.01);
        //printf("%f\n", z[1]);
        x = z[0];
        y = z[1];
        t = i * 1e-2;
        fprintf(fp, "%7.6e %7.6e %7.6e\n", t, x, y);
    }
    fclose(fp);

    /*Solve the Brusselator with adapted Bulirsch-Stoer method*/
    x = 0;
    y = 0;
    t = t_0;
    lev = (int *) malloc(32 * sizeof(int));
    lev[0] = 0;
    lev[1] = 0;
    i = 1;

    sprintf(filename, "BZ_reaction_adapted-BS.txt");
    fp = fopen(filename, "w");
    fprintf(fp, "The Belousov–Zhabotinsky reaction with adapted Bulirsch-Stoer method\n");
    fprintf(fp, "dx / dt = 1 - (b + 1) * x + a * x ^ 2 * y\n");
    fprintf(fp, "dy / dt = b * x - a * x ^ 2 * y\n");
    fprintf(fp, "a = 1, b = 3, x_0 = 0, y_0 = 0, accuracy delta = 1e-10, H = 20\n\n");
    fprintf(fp, "       t            x             y\n");
    fprintf(fp, "%7.6e %7.6e %7.6e\n", t, x, y);

    //z = Bulirsch_Stoer_2D(x, y, t, fx, fy, 20);
    //printf("%e %e a\n", z[0], z[1]);
    //printf("%e err\n", error);

    do
    {
        h_i = H / pow(2, lev[i]);
        z = Bulirsch_Stoer_2D(x, y, t, fx, fy, h_i);
        //printf("%e %e\n", h_i, error);
        if(error >= h_i * delta)
        {
            i++;
            lev[i] = lev[i - 1] + 1;
            lev[i - 1] = lev[i];
        }
        else if(error < h_i * delta)
        {
            t += h_i;
            i--;
            x = z[0];
            y = z[1];
            fprintf(fp, "%7.6e %7.6e %7.6e\n", t, x, y);
            //printf("%e\n", t);
        }
    } while (i > 0);

    fclose(fp);
    free(z);
}

double fx(double x, double y, double t)
{
    return (1 - (b + 1) * x + a * x * x * y);
}

double fy(double x, double y, double t)
{
    return (b * x - a * x * x * y);
}

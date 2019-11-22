/* Hw_1-3.c */

/*
In cosmology, density fluctuations in the matter distribution are characterized by a power spectrum P(k),
the rms (root-mean-square) amplitude fluctuations of the density waves,
as a function of wavenumber k (with units of (h / Mpc)).
In configuration space, these density fluctuations are described by correlation function \xi(r),
at a given scale r, usually in Mpc/h.
These two are quantities related by

\xi(r) = (1 / (2 * \pi ^2)) * Integrate [k ^ 2 * P(k) *(sin(k * r) / k * r), {k, 0, \infinity}]

Attached is a tabulated power spectrum, lcdm_z0.matter_pk .
The first column is k, and the second column is P(k).
Using whatever integration method preferred,
use the equation above to calculate \xi(r) in the range r = [50 , 120] (Mpc / h).
The power spectrum is tabulated in logarithmic intervals in k, due to its power-law like nature.
An interpolation technique may be chosen to help evaluate the integral.

Around k ~ 0.1, you can see oscillatory behavior in P(k). These are called the "baryon wiggles",
and they manifest as a single "bump" in the correlation function at large scales.
Using your calculation for \xi(r), determine the scale r of the peak of this bump.
Make a plot of r ^ 2 * \xi(r) over the required range in r (multiplying by r ^ 2 visually enhances the bump).
Indicate on this plot the scale of the peak, also known as the "baryon acoustic oscillation" (BAO) peak

Notes:
a) You may use a prepackaged routine (or Numerical Recipe code) if you use spline interpolation.
b) Formally, the limits of the integral are from k = 0 to k = \infinity. Note that P(0) = 0.
   You may choose a finite upper limit, provided that you can determine if your limit is robust.
c) the 'h' in the distance units refers to the Hubble constant h = H_0 / 100, which sets the distance scale,
   and is thus incorporated into the distance units, since its value is unknown.
*/

#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <stdbool.h> //standard boolean library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routines

#include "Pk.h"      //P(k) header file

void main()
{
    //return char* will lead to an error in the last digit of the string, instead, return an array char []
    struct cspline cspline_ntr_bc(struct two_dim_data f);    //the second derivatives of natural cubic spline's result, can be used to calculate the interpolation function

    /*Input data, either manually or by file*/
    printf("Calculate correlation function \\ksi(r) from power spectrum P(k)\n\n");
    input_Pk();
    //for(i = 0; i <= N-1; i++) printf("%7.6e\n", Pk[i].P);

    /*Calculate log(k)-log(P) for cubic spline*/
    strcpy(logPk.name, "log_of_");
    strcat(logPk.name, fn.filename);
    //printf("%s", logPk.name);                              //should be a static array instead of a dynamical to operate well
    logPk.n = N - 1;                                         //total number of x_0, ..., x_n is N = n + 1
    logPk.x = (float *) malloc((logPk.n + 1) * sizeof(float));
    logPk.y = (float *) malloc((logPk.n + 1) * sizeof(float));
    for(int i = 0; i <= logPk.n; i++)
        logPk.x[i] = log10(Pk[i].k);
    for(int i = 0; i <= logPk.n; i++)
        logPk.y[i] = log10(Pk[i].P);
    //printf("%7.6e\n", logPk.x[0]);

    /*Calculate the cubic spline interpolation of log(k)- log(P)*/
    cs_logPk = cspline_ntr_bc(logPk);                        //the coefficient matrix of second derivative vector is strictly diagonally dominant, so there should always be an answer
    //printf("%7.6e %7.6e %7.6e %7.6e\n", cs_logPk.a[cs_logPk.n - 1], cs_logPk.b[cs_logPk.n - 1], cs_logPk.c[cs_logPk.n - 1], cs_logPk.y[cs_logPk.n - 1]);

    /*Calculate \xi(r) for different r*/
    corr_func(Pk, cs_logPk);   //calculate \xi(r) in the user defined range of r and store in a file

    free(logPk.x);
    free(logPk.y);
    free(cs_logPk.x);
    free(cs_logPk.y);
    free(cs_logPk.m);
    free(cs_logPk.a);
    free(cs_logPk.b);
    free(cs_logPk.c);
    free(cs_logPk.h);
}

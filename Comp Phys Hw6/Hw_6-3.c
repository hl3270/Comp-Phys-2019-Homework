/*Hw_6-3.c*/

/*
Power spectra of random numbers and random walks:

a) Construct your own linear congruential random number generator to generate uniformly distributed random numbers between [0, 1).
You can choose your values for the constants either from Newman or from the Wikipedia page on LCG.
Make sure you clearly state your choice in your writeup.

b) Using this LCG, create a function that generate Gaussian random variables. Demonstrate that your code works
by comparing a histogram of 10,000 generated values to the unit Gaussian (zero mean, unit standard deviation).
Make the y-axis a log scale - this more easily shows that the frequency of rare events is correctly calculated in your code.

c) Using your code from homework #3, for producing the discrete Fourier transform of a set of data, calculate the power spectrum of your list of 10,000 random Gaussian numbers.
Show that this list has the correct scaling with wavenumber k with a plot showing log(P) vs. log(k).

d) Use this list of random Gaussian numbers to construct a random walk - i.e.,
the i-th value of the sequence is (i - 1)-th value in the sequence plus the i-th Gaussian number. Show a plot of your random walk as a function of iteration i.

e) Calculate the power spectrum of the random walk, and demonstrate that it has the correct scaling with k, using a plot of log(P) vs. log(k).
(Hint: you should have k ^ (- 2)). Remember that the power spectrum is the square of the magnitude of the Fourier coefficients.
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines

#define PI (3.14159265)

void main()
{
    struct complex_array dft(int N, double *y);

    int N;                                //amount of random numbers to generate
    long long int a;                      //LCG : x' = (a * x + c) % m
    long long int c;
    long long int m;
    long long int x;
    long long int x0;                     //initial of x
    double random;                        //uniform random numbers within [0, 1)
    double sigma;
    double r;
    double theta;
    double *g;                            //an array of Gaussian random numbers

    int *hist;                            //histogram array
    int nhist;                            //number of intervals in histogram
    double lhist;                         //left end of the histogram
    double whist;                         //bar width of the histogram

    double *rwalk;                        //random walk

    double *real;                         //real part of the Fourier coefficients
    double *imaginary;                    //imaginary part of the Fourier coefficients
    double *ck2;                          //power spectrum

    int i;
    int j;
    int k;
    int n;
    FILE *fp;
    char filename[64];
    FILE *fp1;
    char filename1[64];

    printf("Power spectra of random numbers and random walks\n\n");

    /*Get N random Gaussian numbers*/
    N = 10000;
    sigma = 1;

    g = (double *) malloc(N * sizeof(double));

    sprintf(filename1, "rand_uni_N_%d.txt", N);
    fp1 = fopen(filename1, "w");
    fprintf(fp1, "Uniform random numbers N = %d\n\n", N);

    sprintf(filename, "rand_Gauss_N_%d.txt", N);
    fp = fopen(filename, "w");
    fprintf(fp, "Random Gaussian numbers N = %d\n\n", N);

    a = 1664525;
    c = 1013904223;
    m = 4294967296;                       //2 ^ 32
    x0 = (long long int) time(NULL);      //choose the starting point as the time
    //x0 = 1578719623;
    //printf("%lld\n", x0);
    x = x0;
    for(i = 0; i < N; i++)               //linear congruential random number generator within [0, 1) and at the same time calculate the random Gaussian number
    {
        x = (a * x + c) % m;              //calculate r
        random = (double) x / m;
        //printf("%lld\n", x);
        //printf("%lf\n", random);
        fprintf(fp1, "%7.6e\n", random);
        r = sqrt(- 2 * sigma * sigma * log(1 - random));
        //printf("%f\n", r);

        x = (a * x + c) % m;              //calculate theta
        random = (double) x / m;
        fprintf(fp1, "%7.6e\n", random);
        theta = 2 * PI * random;

        g[i] = r * cos(theta);            //random Gaussian number array
        //printf("%7.6e\n", g[i]);
        fprintf(fp, "%7.6e\n", g[i]);
    }
    fclose(fp1);
    fclose(fp);

    /*Histogram of random Gaussian numbers*/
    nhist = 80;                           //histogram within [-4, 4), with width 0.1
    lhist = - 4;
    whist = 0.1;

    hist = (int *) malloc(nhist * sizeof(int));
    for(i = 0; i < nhist; i++)
        hist[i] = 0;
    for(i = 0; i < N; i++)
    {
        j = (int) ((g[i] - lhist) / whist);
        hist[j] += 1;
    }
    //for(i = 0; i < nhist; i++)
        //printf("%d\n", hist[i]);
    sprintf(filename, "hist_rand_Gauss_N_%d.txt", N);
    fp = fopen(filename, "w");
    fprintf(fp, "Histogram of random Gaussian numbers\n\n");
    fprintf(fp, "N = %d, range [%3.2e,%3.2e], %d bins, bin width %3.2e\n\n", N, lhist, (lhist + nhist * whist), nhist, whist);
    for(i = 0; i < nhist; i++)
    {
        fprintf(fp, "%d\n", hist[i]);
    }
    fclose(fp);
    free(hist);

    /*Calculate the power spectrum of the random Gaussian numbers*/
    sprintf(filename, "pow_spec_rand_Gauss_N_%d.txt", N);
    fp = fopen(filename, "w");
    fprintf(fp, "Power spectrum of random Gaussian numbers N = %d\n\n", N);
    fprintf(fp, "   k       |c_k|^2\n");

    real = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    imaginary = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    ck2 = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    for(k = 0; k < (N / 2 + 1); k++)
    {
        real[k] = 0;
        imaginary[k] = 0;
        ck2[k] = 0;
    }

    for(k = 0; k < (N / 2 + 1); k++)
    {
        for(n = 0; n <= N - 1; n++)
        {
            real[k] += g[n] * cos(2 * PI * k * n / N);
            imaginary[k] += g[n] * sin(2 * PI * k * n / N);
        }
    }
    for(k = 0; k < (N / 2 + 1); k++)
    {
        ck2[k] = real[k] * real[k] + imaginary[k] * imaginary[k];
        //printf("%f\n", ck2[k]);
        fprintf(fp, "%6d %7.6e\n", k, ck2[k]);
    }
    fclose(fp);
    free(real);
    free(imaginary);
    free(ck2);

    /*Random walk*/
    sprintf(filename, "rand_walk_N_%d.txt", N);
    fp = fopen(filename, "w");
    fprintf(fp, "Random walk by random Gaussian numbers N = %d\n\n", N);

    rwalk = (double *) malloc(N * sizeof(double));
    rwalk[0] = g[0];
    //printf("%7.6e\n", rwalk[0]);
    fprintf(fp, "%7.6e\n", rwalk[0]);
    for(i = 1; i < N; i++)
    {
        rwalk[i] = rwalk[i - 1] + g[i];
        //printf("%7.6e\n", rwalk[i]);
        fprintf(fp, "%7.6e\n", rwalk[i]);
    }
    fclose(fp);

    /*Calculate the power spectrum of the random walk*/
    sprintf(filename, "pow_spec_rand_walk_N_%d.txt", N);
    fp = fopen(filename, "w");
    fprintf(fp, "Power spectrum of random walk N = %d\n\n", N);
    fprintf(fp, "   k       |c_k|^2\n");

    real = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    imaginary = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    ck2 = (double *) malloc((int)(N / 2 + 1) * sizeof(double));
    for(k = 0; k < (N / 2 + 1); k++)
    {
        real[k] = 0;
        imaginary[k] = 0;
        ck2[k] = 0;
    }

    for(k = 0; k < (N / 2 + 1); k++)
    {
        for(n = 0; n <= N - 1; n++)
        {
            real[k] += rwalk[n] * cos(2 * PI * k * n / N);
            imaginary[k] += rwalk[n] * sin(2 * PI * k * n / N);
        }
    }
    for(k = 0; k < (N / 2 + 1); k++)
    {
        ck2[k] = real[k] * real[k] + imaginary[k] * imaginary[k];
        //printf("%f\n", ck2[k]);
        fprintf(fp, "%6d %7.6e\n", k, ck2[k]);
    }
    fclose(fp);
    free(real);
    free(imaginary);
    free(ck2);

    free(rwalk);
    free(g);
}

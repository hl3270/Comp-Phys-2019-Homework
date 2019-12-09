/*Hw_3-1.c*/

/*
Exercise 7.2 (Computational Physics, Mark Newman)
Detecting periodicity

In the on-line resources there is a file called "sunspots.txt" ,which contains the observed number of sunspots
on the Sun for each month since January 1749.
The file contains two columns of numbers, the first representing the month and the second being the sunspot number.

a) Write a program that reads the data in the file and makes a graph of sunspots as a function of time.
You should see that the number of sunspots has fluctuated on a regular cycle for as long as observations have been recorded.
Make an estimate of the length of the cycle in months.

b) Modify your program to calculate the Fourier transform of the sunspot data and then make a graph of the magnitude squared
|c_k|2 of the Fourier coefficients as a function of k - also called the power spectrum of the sunspot signal.
You should see that there is a noticeable peak in the power spectrum at a nonzero value of k.
The appearance of this peak tells us that there is one frequency in the Fourier series that has a higher amplitude
than the others around it - meaning that there is a large sine-wave term with this frequency,
which corresponds to the periodic wave you can see in the original data.

c) Find the approximate value of k to which the peak corresponds. What is the period of the sine wave with this value of k?
You should find that the period corresponds roughly to the length of the cycle that you estimated in part (a).
This kind of Fourier analysis is a sensitive method for detecting periodicity in signals.
Even in cases where it is not clear to the eye that there is a periodic component to a signal,
it may still be possible to find one using a Fourier transform.

You are required to program your own DFT code and use it on the sunspot data.
*/

#include <stdlib.h>        //standard library
#include <stdio.h>         //standard input/output library
#include <string.h>        //string routines
#include <math.h>          //math routines
#include <malloc.h>        //memory allocation routine
#include <stdbool.h>       //standard boolean routine

#include "data_file.h"
#include "fourier_transform.h"

#define PI (3.14159265)

void main()
{
    FILE *fp;
    int N;                                            //number of data points
    int *month;
    float *sunspot;                                   //number of sunspots every month
    struct complex_array dft_sunspot;                 //Fourier coefficients of DFT of sunspots data


    //printf("Discrete Fourier Transform\n\n");
    //printf("Enter a data file for DFT : ");
    //scanf("%s", fn.filename);
    //fp = fopen(fn.filename, "r");
    fp = fopen("sunspots.txt", "r");                  //read data from the file
    N = row_file("sunspots.txt");
    month = (int *) malloc(N * sizeof(int));
    sunspot = (float *) malloc(N * sizeof(float));
    for(int i = 1; i <= N; i++)
    {
        fscanf(fp, "%d %f\n", &month[i], &sunspot[i]);
        //printf("%d %f\n", month[i],sunspot[i]);
    }

    sprintf(data_name, "sunspots");
    dft_sunspot = dft(N, sunspot);

    fp = fopen("dft_mag_sq_sunspots.txt", "w");
    fprintf(fp, "Magnitude square of the discrete Fourier transform of sunspots.txt\n\n");
    fprintf(fp, "  k        |c_k|^2    \n");
    printf("Magnitude square of the discrete Fourier transform of sunspots.txt\n\n");
    printf("  k        |c_k|^2    \n");
    for(int k = 0; k < (N / 2 + 1); k++)              //power spectrum
    {
        printf("%5d %18f\n", k, (dft_sunspot.real[k] * dft_sunspot.real[k] + dft_sunspot.imaginary[k] * dft_sunspot.imaginary[k]));
        fprintf(fp, "%5d %18f\n", k, (dft_sunspot.real[k] * dft_sunspot.real[k] + dft_sunspot.imaginary[k] * dft_sunspot.imaginary[k]));
    }

    fclose(fp);

    free(month);
    free(sunspot);
    free(dft_sunspot.real);
    free(dft_sunspot.imaginary);
}

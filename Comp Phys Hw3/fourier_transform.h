/*
fourier_transform.h

Functions of Fourier transforms : DFT, FFT, ...

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routine
*/

#define PI (3.14159265)

char data_name[64];      //name of data for Fourier transform

struct complex_array     //structure of an N-element complex number array
{
    int N;
    float *real;
    float *imaginary;
};

/*
Discrete Fourier transform

c_k = sum_{n = 0...N - 1} y_n * exp(- i * 2 * PI * k * n / N)

c_{N - r) = c_r*

|c_k| ^ 2 = (sum_{n = 0...N - 1} y_n * cos(2 * PI * k * n / N)) ^ 2 + (sum_{n = 0...N - 1} y_n * cos(2 * PI * k * n / N)) ^ 2
*/
struct complex_array dft(int N, float *y)
{
    struct complex_array dFT;                                               //result coefficients of discrete Fourier transform
    FILE *fp;

    dFT.N = N;                                                              //number of data
    dFT.real = (float *) malloc((int)(N / 2 + 1) * sizeof(float));          //real part of the Fourier coefficient array
    dFT.imaginary = (float *) malloc((int)(N / 2 + 1) * sizeof(float));     //imaginary part of the Fourier coefficient array

    for(int k = 0; k < (N / 2 + 1); k++)
    {
        for(int n = 0; n <= N - 1; n++)
        {
            dFT.real[k] += y[n] * cos(2 * PI * k * n / N);
            dFT.imaginary[k] += y[n] * sin(2 * PI * k * n / N);
        }
    }

    sprintf(fn.filename, "dft_%s.txt", data_name);
    fp = fopen(fn.filename, "w");
    fprintf(fp, "Discrete Fourier Transform of %s\n\n", data_name);
    fprintf(fp, "  k         Re(c_k)          Im(c_k)     \n");
    for(int k = 1; k < (N / 2 + 1); k++)
    {
        fprintf(fp, "%5d %16f %16f\n", k, dFT.real[k], dFT.imaginary[k]);
    }
    fclose(fp);

    return dFT;
}

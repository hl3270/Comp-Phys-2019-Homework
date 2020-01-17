/*Hw_6-1.c*/

/*
Exercise 10.9 (Computational Physics, Mark Newman)
The Ising Model

You do not need to make a movie, but you should show the time evolution of your model in some manner.
A multi-panel plot logarithmically-spaced time intervals is probably most efficient.
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines
#include <time.h>       //time routines

#define SEED time(NULL)

void main()
{
    int s[20][20];        //20 * 20 spin lattice
    int rx;               //choose a random x
    int ry;               //choose a random y
    double r;             //random number
    double *E;            //total energy
    int ss;               //sum of s_i * s_j
    int *M;               //magnetization
    double J;             //interaction constant
    double T;             //temperature
    double kb;            //Boltzmann constant
    int N;                //target Monte Carlo steps
    int n;                //current step
    int i;                //x index
    int j;                //y index
    FILE *fp;
    char filename[64];

    J = 1;
    T = 1;
    kb = 1;
    N = 1e6;

    srand(SEED);

    /*Initialize spin lattice*/
    n = 0;
    for(i = 0; i < 20; i++)
    {
        for(j = 0; j < 20; j++)
            {
                r = (double) rand() / (RAND_MAX + 1);
                //printf("%lf\n", r);
                if(r < 0.5)
                    s[i][j] = - 1;
                else
                    s[i][j] = 1;
                //printf("%d\n", s[i][j]);
            }
    }
    sprintf(filename, "Ising_spin_T_%2.1f_n_%d.txt", T, n);
    fp = fopen(filename, "w");
    fprintf(fp, "Spin lattice of Ising model with T = %2.1f at step n = %d\n\n", T, n);
    for(i = 0; i < 20; i++)
    {
        for(j = 0; j < 20; j++)
        {
            fprintf(fp, "%d ", s[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    /*Initialization of energy and magnetization*/
    E = (double *) malloc((N + 1) * sizeof(double));
    M = (double *) malloc((N + 1) * sizeof(double));
    E[0] = 0;
    M[0] = 0;
    for(i = 0; i < 20; i++)
    {
        for(j = 0; j < 20; j++)
        {
            M[0] += s[i][j];
        }
    }
    //printf("%d\n", M[0]);
    ss = 0;
    for(i = 0; i < 19; i++)
    {
        for(j = 0; j < 19; j++)
        {
            ss += s[i][j] * s[i + 1][j];
            ss += s[i][j] * s[i][j + 1];
        }
    }
    for(i = 0; i < 19; i++)
        ss += s[i][19] * s[i + 1][19];
    for(j = 0; j < 19; j++)
        ss += s[19][j] * s[19][j + 1];
    //printf("%d\n", ss);
    E[0] = - J * ss;
    //printf("%lf\n", E[0]);

    /*Flip spin until target step*/
    for(n = 1; n <= N; n++)
    {
        /*Flip s random spin*/
        rx = (int) (20 * rand() / (RAND_MAX + 1));
        ry = (int) (20 * rand() / (RAND_MAX + 1));
        //printf("%d\n", rx);
        s[rx][ry] = - 1 * s[rx][ry];
        ss = 0;
        for(i = 0; i < 19; i++)
        {
            for(j = 0; j < 19; j++)
            {
                ss += s[i][j] * s[i + 1][j];
                ss += s[i][j] * s[i][j + 1];
            }
        }
        for(i = 0; i < 19; i++)
            ss += s[i][19] * s[i + 1][19];
        for(j = 0; j < 19; j++)
            ss += s[19][j] * s[19][j + 1];
        E[n] = - J * ss;
        //printf("%lf ", E[n]);

        if(E[n] > E[n - 1])
        {
            r = (double) rand() / (RAND_MAX + 1);
            if(r >= exp((E[n - 1] - E[n]) / (kb * T)))  //not accept move
            {
                E[n] = E[n - 1];
                M[n] = M[n - 1];
                s[rx][ry] = - 1 * s[rx][ry];
            }
            else                                        //Boltzmann accept move
            {
                M[n] = 0;
                for(i = 0; i < 20; i++)
                {
                    for(j = 0; j < 20; j++)
                    {
                        M[n] += s[i][j];
                    }
                }
            }
        }
        else                                            //lower energy accept move
        {
            M[n] = 0;
            for(i = 0; i < 20; i++)
            {
                for(j = 0; j < 20; j++)
                {
                    M[n] += s[i][j];
                }
            }
        }
        //printf("%d\n", M[n]);

        if((n == 1e2) || (n == 1e4) || (n == 1e6))
        {
            sprintf(filename, "Ising_spin_T_%2.1f_n_%d.txt", T, n);
            fp = fopen(filename, "w");
            fprintf(fp, "Spin lattice of Ising model with T = %2.1f at step n = %d\n\n", T, n);
            for(i = 0; i < 20; i++)
            {
                for(j = 0; j < 20; j++)
                {
                    fprintf(fp, "%d ", s[i][j]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }
    }

    sprintf(filename, "Ising_mag_T_%2.1f.txt", T);
    fp = fopen(filename, "w");
    fprintf(fp, "Magnetization of Ising model using Markov chain Monte Carlo simulation\n");
    fprintf(fp, "T = %2.1f, step n = [0, %d]\n\n", T, N);
    for(n = 0; n <= N; n++)
    {
        fprintf(fp, "%d\n", M[n]);
    }
    fclose(fp);

    free(E);
    free(M);
}

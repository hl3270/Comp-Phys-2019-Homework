/*Hw_4-2cd.c*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines

int dim = 4;
double r_s = 1e-7;

void main()
{
    double *f(int d, double *x, double t, double A, double B);
    double tSBHB(int d, double *r, double *(*func_ptr)(int d, double *r, double t, double A, double B), double t_0, double delta);

    double *r;               //r = (x, y, v_x, v_y)
    double t_0;
    double H;
    double delta;
    double t;
    double A;
    double B;
    int exit;

    r = (double *) malloc(dim * sizeof(double));
    r[0] = 1;
    r[1] = 0;
    r[2] = 0;
    r[3] = 0.4;
    t_0 = 0;
    delta = 1e-8;

    for(A = 0.5; A <= 10; A += 0.5)
    {
        for(j = 1; j <= 2; B += 0.5)
        {
            A = i / 2;
            B = j / 2;
            r[0] = 1;
            r[1] = 0;
            r[2] = 0;
            r[3] = 0.4;
            t = tSBHB(dim, r, f, t_0, delta);
            printf("%e %e %e %e\n", A, B, (B / A), t);
        }
    }
}

/*
Calculate the time it takes for SBHB to reach Schwarzschild radius
*/
double tSBHB(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t_0, double delta)
{
    double *RK4_step(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t, double h);

    double *r_1;
    double *r_2;
    double t;
    double h;
    double delx;
    double rho;
    double r_peri;
    int i;

    h = 1e-3;
    r_1 = (double *) malloc(d * sizeof(double));
    r_2 = (double *) malloc(d * sizeof(double));

    do
    {
        do
        {
            r_2 = RK4_step(d, r, func_ptr, t, h);
            r_1 = RK4_step(d, r_2, func_ptr, t + h, h);
            r_2 = RK4_step(d, r, func_ptr, t, 2 * h);

            delx = 0;
            for(i = 0; i < 2; i++)   //??
            {
                delx  = delx + (r_2[i] - r_1[i]) * (r_2[i] - r_1[i]);
            }
            delx = sqrt(delx);
            rho = 30 * h * delta / delx;

            if(rho > 1)
            {
                t += 2 * h;
                printf("%e tt\n", t);
                for(i = 0; i < d; i++)
                {
                    r[i] = r_1[i];
                }
            }

            h = h * sqrt(sqrt(rho));
        } while(rho <= 1);

        r_peri = sqrt(r[0] * r[0] + r[1] * r[1]);
    } while(r_peri > r_s);

    return t;

    free(r_1);
    free(r_2);
}

double *RK4_step(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t, double h)
{
    double *k_1;
    double *k_2;
    double *k_3;
    double *k_4;
    double *r_temp;
    int i;

    k_1 = (double *) malloc(d * sizeof(double));
    k_2 = (double *) malloc(d * sizeof(double));
    k_3 = (double *) malloc(d * sizeof(double));
    k_4 = (double *) malloc(d * sizeof(double));
    r_temp = (double *) malloc(d * sizeof(double));

    k_1 = (*func_ptr)(d, r, t);
    for(i = 0; i < d; i++)
    {
        k_1[i] = h * k_1[i];
        r_temp[i] = r[i] + k_1[i] / 2;
    }
    k_2 = (*func_ptr)(d, r_temp, t + h / 2);
    for(i = 0; i < d; i++)
    {
        k_2[i] = h * k_2[i];
        r_temp[i] = r[i] + k_2[i] / 2;
    }
    k_3 = (*func_ptr)(d, r_temp, t + h / 2);
    for(i = 0; i < d; i++)
    {
        k_3[i] = h * k_3[i];
        r_temp[i] = r[i] + k_3[i];
    }
    k_4 = (*func_ptr)(d, r_temp, t + h);
    for(i = 0; i < d; i++)
    {
        k_4[i] = h * k_4[i];
    }

    for(i = 0; i < d; i++)
    {
        r_temp[i] = r[i] + (k_1[i] + 2 * k_2[i] + 2 * k_3[i] + k_4[i]) / 6;
    }

    return r_temp;

    free(k_1);
    free(k_2);
    free(k_3);
    free(k_4);
}

double *f(int d, double *x, double t)
{
    double *y;
    double r;
    double v;

    y = (double *) malloc(d * sizeof(double));

    r = sqrt(x[0] * x[0] + x[1] * x[1]);
    v = sqrt(x[2] * x[2] + x[3] * x[3]);

    y[0] = x[2];
    y[1] = x[3];
    y[2] = - x[0] / (4 * r * r * r) - A * x[2] / (B + v * v * v);
    y[3] = - x[1] / (4 * r * r * r) - A * x[3] / (B + v * v * v);

    return y;
}

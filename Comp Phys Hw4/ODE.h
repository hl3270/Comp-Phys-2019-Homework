/*
ODE.h

Methods of solving ordinary differential equations.

(Here to be simple, only list codes for two dimension, which is used in the homeworks. Generalization is to be followed)

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <stdbool.h> //standard boolean library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routines
*/

#include "SBHB.h"

int N = 8;              //maximum number for modified midpoint method
double H = 20;
double t_0 = 0;
double error;
double delta = 1e-10;

/*
Modified midpoint method for two dimension

Calculate x(t+H0) from x(t) with n steps
*/
double *mod_mid_2D(double x, double y, double t, double (*fx_ptr)(double x, double y, double t), double (*fy_ptr)(double x, double y, double t), int n, double H0)
{
    int m;
    double x1n;
    double y1n;
    double x2n;
    double y2n;
    double h;

    h = H0 / n;
    x2n = x + h * fx_ptr(x, y, t);
    y2n = y + h * fy_ptr(x, y, t);
    x1n = x + h * fx_ptr(x, y, t + h / 2);
    y1n = y + h * fy_ptr(x, y, t + h / 2);

    for(m = 1; m <= n - 1; m++)
    {
        x2n = x2n + h * fx_ptr(x1n, y1n, (t + m * h));
        y2n = y2n + h * fy_ptr(x1n, y1n, (t + m * h));
        x1n = x1n + h * fx_ptr(x2n, y2n, (t + (m + 1 / 2) * h));
        y1n = y1n + h * fy_ptr(x2n, y2n, (t + (m + 1 / 2) * h));
        //printf("%f %d\n", x2n, m);
    }

    x2n = (x1n + x2n + h * fx_ptr(x1n, y1n, (t + H0)) / 2) / 2;
    y2n = (y1n + y2n + h * fy_ptr(x1n, y1n, (t + H0)) / 2) / 2;

    double *z;
    z = (double *) malloc(2 * sizeof(double));
    z[0] = x2n;
    z[1] = y2n;
    //printf("%f\n", z[0]);

    return z;
}

/*
Bulirsch-Stoer Method for two dimension
*/
double *Bulirsch_Stoer_2D(double x, double y, double t, double (*fx_ptr)(double x, double y, double t), double (*fy_ptr)(double x, double y, double t), double H0)
{
    int n;
    int m;
    int i;
    float j;
    double *rx;  //storage of Richardson extrapolation
    double *ry;
    double *r;
    double ex;
    double ey;
    double *mod_mid_2D(double x, double y, double t, double (*fx_ptr)(double x, double y, double t), double (*fy_ptr)(double x, double y, double t), int n, double H0);

    error = 2 * H * delta;
    r = (double *) malloc(2 * sizeof(double));
    rx = (double *) malloc(N * (N + 1) * sizeof(double) / 2);
    ry = (double *) malloc(N * (N + 1) * sizeof(double) / 2);

    /*Initialization*/
    r = mod_mid_2D(x, y, t, fx_ptr, fy_ptr, 1, H0);
    rx[0] = r[0];
    ry[0] = r[1];

    /*Extrapolation*/
    for(n = 2; n <= N; n++)
    {
        for(m = 1; m <= n; m++)
        {
            if(m == 1)
            {
                r = mod_mid_2D(x, y, t, fx_ptr, fy_ptr, n, H0);
                rx[n * (n - 1) / 2] = r[0];
                ry[n * (n - 1) / 2] = r[1];
                //printf("%e %e i\n", rx[n * (n - 1) / 2], ry[n * (n - 1) / 2]);
            }

            else
            {
                j = pow(n, (2 * m - 2)) / pow((n - 1), (2 * m - 2)) - 1;
                //printf("%f j \n", j);
                rx[n * (n - 1) / 2 + m - 1] = rx[n * (n - 1) / 2 + m - 2] + (rx[n * (n - 1) / 2 + m - 2] - rx[(n - 1) * (n - 2) / 2 + m - 2]) / j;
                ry[n * (n - 1) / 2 + m - 1] = ry[n * (n - 1) / 2 + m - 2] + (ry[n * (n - 1) / 2 + m - 2] - ry[(n - 1) * (n - 2) / 2 + m - 2]) / j;
                //printf("%e %e ii\n", rx[n * (n - 1) / 2 + m - 1], ry[n * (n - 1) / 2 + m - 1]);
            }
        }

        ex = (rx[n * (n + 1) / 2 - 2] - rx[n * (n - 1) / 2 - 1]) / j;
        ey = (ry[n * (n + 1) / 2 - 2] - ry[n * (n - 1) / 2 - 1]) / j;
        //printf("%e %e exey\n", ex, ey);
        error = sqrt(ex * ex + ey * ey);
        //printf("%e e\n", error);

        if(error >= 1e30)
        {
            r[0] = rx[0];
            r[1] = ry[0];
            error = 2 * H0 * delta;
            break;
        }
        else if(error < H0 * delta)
        {
            r[0] = rx[n * (n + 1) / 2 - 1];
            r[1] = ry[n * (n + 1) / 2 - 1];
            break;
        }
    }

    //printf("%e r\n", r[0]);
    return r;

    free(rx);
    free(ry);
}

/*
4th order Runge-Kutta method for d dimension of one step h

k_1 = h * f (r, t)
k_2 = h * f(r + k_1 / 2, t + h / 2)
k_3 = h * f(r + k_2 / 2, t + h / 2)
k_4 = h * f(r + k_3, t + h)
r (t + h) = r (t) + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6

Here r is d dimensional
*/
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

    //printf("%e %e %e %e\n", k_1[0], k_1[1], k_1[2], k_1[3]);
    //printf("%e %e %e %e\n", k_2[0], k_2[1], k_2[2], k_2[3]);
    //printf("%e %e %e %e\n", k_3[0], k_3[1], k_3[2], k_3[3]);
    //printf("%e %e %e %e\n", k_4[0], k_4[1], k_4[2], k_4[3]);

    return r_temp;

    free(k_1);
    free(k_2);
    free(k_3);
    free(k_4);
}

/*
4th order Runge-Kutta method for d dimension
with step size h, from t = t_0 to t = t_0 + H
*/


/*
4th order Runge-Kutta method for d dimension with adaptive time step
*/
void adapt_RK4(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t_0, double H, double delta)
{
    double *RK4_step(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t, double h);

    FILE *fp;
    char filename[64];
    double *r_1;
    double *r_2;
    double t;
    double h;
    double delx;
    double rho;
    double r_peri;
    int i;

    sprintf(filename, "SBHB_v_%4.3e_A_%3.1f_B_%3.1f_delta_%4.3e.txt", r[3], A, B, delta);
    fp = fopen(filename, "w");
    fprintf(fp, "Solution of SBHB orbit\n\n");
    fprintf(fp, "v_0 = %4e, A = %3lf, B = %3lf, delta = %4e\n\n", r[3], A, B, delta);
    fprintf(fp, "       t           rx             ry            vx             vy             2h\n");

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
            //printf("%e\n", delx);
            rho = 30 * h * delta / delx;
            //printf("%e rho\n", rho);

            if(rho > 1)
            {
                t += 2 * h;
                fprintf(fp, "%7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n", t, r_1[0], r_1[1], r_1[2], r_1[3], 2 * h);

                for(i = 0; i < d; i++)
                {
                    r[i] = r_1[i];
                }
            }

            h = h * sqrt(sqrt(rho));
        } while(rho <= 1);

        r_peri = sqrt(r[0] * r[0] + r[1] * r[1]);
        if(r_peri <= r_s)
        {
            break;
        }

    } while(t <= H);

    fclose(fp);

    free(r_1);
    free(r_2);
}

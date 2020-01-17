/*Hw_4-2.c*/

/*
As per the class notes, you will integrate the orbit of a supermassive black hole binary
in the presence of dynamical friction. The equation of the motion for a BH is :

d^2(\vec{r_BH}) / dt^2 = -G * M_BH / (4 * r_BH ^ 3) * \vec(r_BH) + \dot{\vec{v_DF}}

where the last term on the right hand side is the force from dynamical friction (DF).
We use the following approximation formula for this force :

\dot{\vec{v_DF}} = - A / (v_BH ^ 3 + B) * \vec{v_DF}

where the two constants A and B represent the product of the BH mass and the stellar density
and the velocity dispersion (cubed) of the stellar field, respectively.

The goal of this problem is to determine whether the BH binary can get close enough in pericenter
such that it can lose energy via gravitational radiation.
This happens when r_peri \sim the Schwartschild radius.
Set M = G = 1 for the calculation, and the unit of distance to be 100 parsecs. In these units, r_s is 10^(-7).

You will use your fourth-order Runge-Kutta integrator, with adaptive time steps, to solve this calculation.
Since the orbit is symmetric (i.e. the BHs are the same), you only need to solve the position of one BH.
For all calculations with dynamical friction, make the initial position of the BH x = 1, y = 0,
and the initial velocity equal to 0.8 of that of a circular orbit.

a) First, you need to set your value of d, the error tolerance per unit of time.
To do this, first solve for the orbit without dynamical friction, setting the initial velocity
such that r_peri is 10^(-7), and demonstrating that there is no appreciable loss of accuracy over at least 10 orbits.

b) Solve for the BH orbit with A = B = 1. Make a plot showing the path of the BH in its orbit,
and another showing the log(r) as a function of time. You can stop the integration when r = 10^(-7).

c) Determine how the time it takes to reach the Schwartzchild radius depends on the ratio of B / A.
Discuss (qualitatively) why it looks this way. Include a plot showing this, and convert the time axis to Myr.
Take as your range of A and B values [0.5, 10].

d) Do your results depend on the value of the initial velocity? Do they only depend on the ratio of B / A,
or on the individual values of B and A?
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines

#include "ODE.h"

void main()
{
    double *f(int d, double *x, double t);
    double *RK4_step(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t, double h);
    void adapt_RK4(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t_0, double H, double delta);
    double tSBHB(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t_0, double delta);

    double *r;               //r = (x, y, v_x, v_y)
    double t_0;              //initial time
    double H;                //time range
    double delta;            //error tolerance
    double t;                //present time
    int i = 1;
    int j = 1;
    FILE *fp;
    char filename[64];
    int exit;

    printf("Integration of SBHB orbit using adaptive time step RK4 method\n\n");

    /*sprintf(filename, "t_rs_v_%4.3e.txt", 0.4);
    fp = fopen(filename, "w");
    fprintf(fp, "Time for SBHB to reach Schwartzschild radius for different A and B\n");
    fprintf(fp, "v = %7.6e delta = %7.6e\n\n", 0.4, 1e-8);
    fprintf(fp, "      A              B            B / A         t\n");*/

    beginning:
        r = (double *) malloc(dim * sizeof(double));
        //A = 0.5 * i;
        //B = 0.5 * j;
        A = 0;
        B = 0;
        r[0] = 1;
        r[1] = 0;
        r[2] = 0;
        //r[3] = sqrt(1 / (2 * (1e7 + 1)));
        r[3] = 0.4;
        t_0 = 0;
        H = 50;
        delta = 1e-8;

        //x_temp = f(dim, x, t);
        //printf("%f %f %f %f\n", x_temp[0], x_temp[1], x_temp[2], x_temp[3]);
        //x_temp = RK4_step(dim, x, f, t, 1e-2);
        //printf("%f %f %f %f\n", x_temp[0], x_temp[1], x_temp[2], x_temp[3]);

        /*Calculate and save the orbit for some A, B, and v_0*/
        adapt_RK4(dim, r, f, t_0, H, delta);

        /*Calculate the time for SBHB to reach Schwartzschild radius for different A's and B's*/
        //t = tSBHB(dim, r, f, t_0, delta);
        //fprintf(fp, "%7.6e %7.6e %7.6e %7.6e\n", A, B, (B / A), t);
        //printf("%e %e %e %e\n", A, B, (B / A), t);

        /*i++;
        if(i <= 20)
            goto beginning;

        if(j <= 20)
        {
            i = 1;
            j++;
            goto beginning;
        }*/

        //fclose(fp);

        /*printf("\nTo continue, type 1; \nTo exit, press any other key.\n");  //choose to continue or exit
        scanf("%d", &exit);
        if (exit==1)
        {
            printf("\n");
            goto beginning;
        }
        else exit;*/
}

/*
Calculate the time for SBHB to reach Schwartzschild radius
*/
double tSBHB(int d, double *r_0, double *(*func_ptr)(int d, double *r, double t), double t_0, double delta)
{
    double *RK4_step(int d, double *r, double *(*func_ptr)(int d, double *r, double t), double t, double h);

    double *r;
    double *r_1;
    double *r_2;
    double t;
    double h;
    double delx;
    double rho;
    double r_peri;
    int i;

    h = 1e-3;
    r = (double *) malloc(d * sizeof(double));
    r_1 = (double *) malloc(d * sizeof(double));
    r_2 = (double *) malloc(d * sizeof(double));

    for(int i = 0; i < d; i++)
        r[i] = r_0[i];

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
                for(i = 0; i < d; i++)
                {
                    r[i] = r_1[i];
                }
            }

            h = h * sqrt(sqrt(rho));
        } while(rho <= 1);

        r_peri = sqrt(r[0] * r[0] + r[1] * r[1]);
    } while(r_peri > r_s);

    //printf("%7.6e t\n", t);
    return t;

    free(r);
    free(r_1);
    free(r_2);
}

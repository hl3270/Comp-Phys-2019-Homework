/*Hw_5-2.c*/

/*
You will use the relaxation method to solve for the potential of a series of charged particles.

In the github under HOMEWORK5 there is a data file names "particles.dat".
Inside this data file is a set of (x, y) coordinates for a series of charged particles, all with charge equal to the electron charge.
The range of coordinates is [0,100]. As with example 9.2 in Newman, these particles are placed in a box that is grounded on all 4 sides.

(a) Using the cloud-in-cell technique, assign the charges to a two-dimensional grid of size [0, M] per side, with M = 100.
The centers of each grid cell (i,j) is (i+0.5,j+0.5). Produce an image of the charge density field.

(b) Use standard relaxation method to solve for Poisson's equation. Produce an image showing the resulting potential field.
Note the number of iterations it takes to converge. Use the same convergence criterion as the example in Newman,
such that the maximum difference for any cell in the grid between the current and prior step is 1.0E-10.

(c) Now use the Gauss-Seidel overrelaxation method to solve for Poisson's equation.
Determine the optimal value of the overrelaxation parameter \omega, using one of the techniques discussed in class
and in the chapter on nonlinear equations. Golden ratio search is a good example.
Find the optimal value of \omega to a precision of 0.001.
Produce a plot showing how your answer for \omega evolves with each step in your minimization process.

Just a reminder that the equation you are solving for is:
\nabla ^ 2 (\phi) = - \rho / \epsilon_0
*/

#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <stdbool.h> //standard boolean library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routines

#include "iteration.h"  //iteration methods header file

#define EPS_0 (8.854187817e-12)  //vacuum permittivity : \epsilon_0 / (F / m)
#define E (-1.602176487e-19)     //electron charge : e / C
#define RESULT_NAME "potential"  //name of the iteration result

void main()
{
    struct two_dim_grid_data relaxation_2D(float coef_plus_x, float coef_minus_x, float coef_plus_y, float coef_minus_y, float (*func_ptr)(int x_i, int y_j, float hx, float hy));
    float func_x_y_hx_hy(int x_i, int y_j, float hx, float hy);

    float (*func_ptr)(int x_i, int y_j, float hx, float hy);
    float hx;
    float hy;
    float coef_minus_x;
    float coef_plus_x;
    float coef_minus_y;
    float coef_plus_y;
    struct two_dim_grid_data potential;     //iteration result

    beginning:
        strcpy(result_name, RESULT_NAME);
        printf("Calculate two dimensional electric potential\n\n");
        n_xy = cloud_in_cell_two_d();

        hx = (n_xy.x[n_xy.n_x] - n_xy.x[0]) / n_xy.n_x;  //grid cell width
        hy = (n_xy.y[n_xy.n_y] - n_xy.y[0]) / n_xy.n_y;
        coef_minus_x = hy * hy / (2 *(hx * hx + hy * hy));
        coef_plus_x = hy * hy / (2 *(hx * hx + hy * hy));
        coef_minus_y = hx * hx / (2 *(hx * hx + hy * hy));
        coef_plus_y = hx * hx / (2 *(hx * hx + hy * hy));
        func_ptr = func_x_y_hx_hy;
        //potential = relaxation_2D(coef_plus_x, coef_minus_x, coef_plus_y, coef_minus_y, func_ptr);  //input conditions, parameters, and calculate the relaxation result
        //printf("\n\n");

        //for(int i = 0; i <= potential.n_x; i++)
        //{
            //for(int j = 0; j <= potential.n_y; j++)
            //{
                //printf("%9.8e %9.8e %9.8e\n", potential.x[i], potential.y[j], potential.z[i][j]);  //save the final results in data file
            //}
        //}

        potential = GS_overrelaxation_2D(coef_plus_x, coef_minus_x, coef_plus_y, coef_minus_y, func_ptr);  //input conditions, parameters, and calculate the relaxation result
        printf("\n\n");

        for(int i = 0; i <= potential.n_x; i++)
        {
            for(int j = 0; j <= potential.n_y; j++)
            {
                //printf("%9.8e %9.8e %9.8e\n", potential.x[i], potential.y[j], potential.z[i][j]);  //save the final results in data file
            }
        }

        free(potential.x);
        free(potential.y);
        free(potential.z);

        printf("\nTo continue, type 1; \nTo exit, press any other key.\n");  //choose to continue or exit
        scanf("%d", &exit);
        if (exit == 1)
        {
            printf("\n");
            goto beginning;
        }
        else exit;

}

float func_x_y_hx_hy(int x_i, int y_j, float hx, float hy)
{
    return (E * n_xy.z[x_i][y_j] / (2 * EPS_0 * (1 / pow(hx, 2) + 1 / pow(hy, 2))));
}

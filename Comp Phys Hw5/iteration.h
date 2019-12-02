/*
iteration.h

Iteration methods for one or multi- dimensions.

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <stdbool.h> //standard boolean library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routines
*/

//#include "multi_dim_data.h"
//#include "data_file.h"
#include "particle_density.h"

//variables
char result_name[64];  //name string for the result of iteration

struct one_dim_Diri_bound_cond  //structure of one dimensional Dirichlet boundary conditions
{
    float x_left;   //left boundary position
    float x_right;  //right boundary position
    float f_left;   //left boundary value
    float f_right;  //right boundary value
};
struct one_dim_Diri_bound_cond one_D_Dir_BC;

struct two_dim_Diri_bound_cond       ///////////
{

};

int one_D_n_cell;   //the number of one dimensional grid cell number

struct multi_dim_grid
{
    int dim;        //dimension
    int *n_cell;    //number of cells of each dimension
};
struct multi_dim_grid nD_grid;

float max_dif_tar;  //target maximum difference between the current and the prior steps of all grid points

//functions
/*
1D relaxation method of boundary value problems

First, replace the first and the second derivative in the ODE with their finite-difference approximation :

f'(x) \approx (f(x + h) - f(x - h)) / (2 * h)
f"(x) \approx (f(x + h) - 2 * f(x) + f(x - h)) / h ^ 2

Second, from the new linear equation derive the related relaxation-method equation in the form :

f*(x) = F(f(x - h), f(x + h), g(x, h))

Here, h is the size of grid cell; f(x) is the function value at x of the present iteration, and f*(x) is the value of
the next iteration; g(x, h) is a constant function of position x and cell size h;
and F is a linear function of f(x - h), f(x + h), and g(h).

Given boundary values, if the spectrum radius of the iteration matrix \rho < 1, then f(x) will converge to the correct answer.

In this function, knowing the coefficients and constant term of the relaxation function

f_star(x) = coef_plus * f(x + h) + coef_minus * f(x - h) + g(x, h)

by inputing the boundary conditions : positions x_left, x_right, related function values f_left, f_right,
the number of grid cells n_cell, and the target maximum difference for any cell between two steps max_dif_tar,
this function will return the solution of function in 2D data type (*x, *y) (as well as the iteration process) saved in a data file.
*/
struct two_dim_data relaxation_1D(float coef_plus, float coef_minus, float (*func_ptr)(float x, float h))
{
    void input_one_D_Dir_BC();  //input boundary conditions with different positions
    void input_one_dim_grid();  //input the number of one dimensional grid cells
    void input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    float h;       //size of grid cells
    int i;         //index of points (x_i, f(x_i))
    float *x;      //position x_i
    float *f;      //function values of the prior step
    float *f_star; //function values of the current step
    float *f_temp; //template of f for swap
    float *dif;    //difference at each grid point between two steps
    float max_dif; //the maximum difference of function values between the current and the prior step of all grid points

    struct two_dim_data relax_1D;  //output result of relaxation method for the 1D boundary value problem
    FILE *fp;                      //output data file for the whole iteration process and the result
    char filename[64];             //filename string

    /*Input parameters of relaxation method of 1D boundary value problems*/  //if they are set at main function, then the three steps could be omitted
    input_one_D_Dir_BC();  //input boundary conditions with different positions
    input_one_dim_grid();  //input the number of one dimensional grid cells
    input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    /*Create a data file saving the iteration (process and) result*/
    relax_1D.n = one_D_n_cell + 1;       //grid point number is one more than cell number for 1D
    strcpy(relax_1D.name, result_name);

    sprintf(filename, "relax_1D_%s.txt", relax_1D.name);
    fp = fopen(filename, "w");
    fprintf(fp, "Iteration result of %s using 1D relaxation method\n\n", relax_1D.name);
    fprintf(fp, "Parameters : \n");
    fprintf(fp, "x_left = %9.8e, f_left = %9.8e\n", one_D_Dir_BC.x_left, one_D_Dir_BC.f_left);
    fprintf(fp, "x_right = %9.8e, f_right = %9.8e\n", one_D_Dir_BC.x_right, one_D_Dir_BC.f_right);
    fprintf(fp, "n_cell = %d, max_dif_tar = %9.8e\n\n", one_D_n_cell, max_dif_tar);

    /*Initialize the function values of all grid points*/
    h = (one_D_Dir_BC.x_right - one_D_Dir_BC.x_left) / one_D_n_cell;  //grid cell width, used in constant terms of relaxation equation

    x = (float *) malloc((one_D_n_cell + 1) * sizeof(float));
    f = (float *) malloc((one_D_n_cell + 1) * sizeof(float));
    f_star = (float *) malloc((one_D_n_cell + 1) * sizeof(float));
    f_temp = (float *) malloc((one_D_n_cell + 1) * sizeof(float));
    dif = (float *) malloc((one_D_n_cell + 1) * sizeof(float));

    for(i = 0; i <= one_D_n_cell; i++)            //grid point positions
        x[i] = i * one_D_Dir_BC.x_right / one_D_n_cell + (one_D_n_cell - i) * one_D_Dir_BC.x_left / one_D_n_cell;
    f[0] = one_D_Dir_BC.f_left;                   //boundary conditions
    f[one_D_n_cell] = one_D_Dir_BC.f_right;
    for(i = 1; i <= one_D_n_cell - 1; i++)        //starting points all set to be 0
        f[i] = 0;

    /*Iterate until reach the target maximum difference*/
    do
    {
        max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

        for(i = 1; i <= one_D_n_cell - 1; i++)    //leave boundary unchanged
        {
            f_star[i] = coef_minus * f[i - 1] + coef_plus * f[i + 1] + (*func_ptr)(x[i], h);
            dif[i] = f_star[i] - f[i];
        }

        for(i = 1; i <= one_D_n_cell - 1; i++)    //get the maximum difference of these grid points for the present and previous step
        {
            if(max_dif < fabs(dif[i]))
                max_dif = fabs(dif[i]);
        }

        for(i = 0; i <= one_D_n_cell; i++)        //renew the f(x) with the new f*(x), swap these two arrays
        {
            f_temp[i] = f_star[i];
            f_star[i] = f[i];
            f[i] = f_temp[i];                     //intermediate process
            //f[i] = f_star[i];                     //direct renew
        }
        //printf("%f\n", f[2]);

        //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);       //add if want to know the intermediate process
        //fprintf(fp, "       x           f(x)       \n");
        //for(i = 0; i <= one_D_n_cell; i++)
            //fprintf(fp, "%9.8e, %9.8e\n", x[i], f[i]);
        //fprintf(fp, "\n");
    } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference

    //printf("%f", max_dif);

    /*Return and save the iteration results*/
    relax_1D.x = (float *) malloc((relax_1D.n) * sizeof(float));      //results  to return
    relax_1D.y = (float *) malloc((relax_1D.n) * sizeof(float));
    for(i = 0; i <= one_D_n_cell; i++)
    {
        relax_1D.x[i] = x[i];
        relax_1D.y[i] = f[i];
    }

    //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);
    fprintf(fp, "       x           f(x)       \n");
    for(i = 0; i <= one_D_n_cell; i++)
        fprintf(fp, "%9.8e %9.8e\n", relax_1D.x[i], relax_1D.y[i]);  //save the final results in data file

    fclose(fp);

    return relax_1D;

    free(x);
    free(f);
    free(f_star);
    free(f_temp);
    free(dif);
};

/*
Two dimensional relaxation method

Similar to one dimension

f_star(x, y) = coef_plus_x * f(x + hx, y) + coef_minus_x * f(x - hx, y) + coef_plus_y * f(x, y + hy) + coef_minus_x * f(x, y - hy) + g(x, y, a, b)

since the step of x and y may not necessarily be the same.
*/
struct two_dim_grid_data relaxation_2D(float coef_plus_x, float coef_minus_x, float coef_plus_y, float coef_minus_y, float (*func_ptr)(int x_i, int y_j, float hx, float hy))
{
    //void input_two_D_Dir_BC();  //input boundary conditions with different positions
    //void input_two_dim_grid();  //input the number of two dimensional grid cells
    void input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    float hx;       //x width of grid cells
    float hy;       //y width of grid cells
    int i;         //index i of points (x_i, y_j, f(x_i, y_j))
    int j;         //index j of points (x_i, y_j, f(x_i, y_j))
    float *x;      //position x_i
    float *y;      //position y_j
    float **f;      //function values of the prior step
    float **f_star; //function values of the current step
    float **f_temp; //template of f for swap    //?
    float **dif;    //difference at each grid point between two steps
    float max_dif; //the maximum difference of function values between the current and the prior step of all grid points

    struct two_dim_grid_data relax_2D;  //output result of relaxation method for the 2D boundary value problem
    FILE *fp;                      //output data file for the whole iteration process and the result
    char filename[64];             //filename string

    int iter;          //number of iterations

    /*Input parameters of relaxation method of 2D boundary value problems*/
    //input_two_D_Dir_BC();  //input boundary conditions with different positions   //how
    //input_one_dim_grid();  //input the number of one dimensional grid cells
    input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    /*Create a data file saving the iteration (process and) result*/
    relax_2D.n_x = n_xy.n_x;
    relax_2D.n_y = n_xy.n_y;
    relax_2D.x = n_xy.x;
    relax_2D.y = n_xy.y;
    strcpy(relax_2D.name, result_name);

    sprintf(filename, "relax_2D_%s.txt", relax_2D.name);
    fp = fopen(filename, "w");
    fprintf(fp, "Iteration result of %s using 2D relaxation method\n\n", relax_2D.name);
    fprintf(fp, "Parameters : \n");
    fprintf(fp, "x \\in [%9.8e, %9.8e], y \\in [%9.8e, %9.8e]\n", relax_2D.x[0], relax_2D.x[relax_2D.n_x], relax_2D.y[0], relax_2D.y[relax_2D.n_y]);
    fprintf(fp, "Boundary value : all 0 \n");  //
    fprintf(fp, "n_x = %d, n_y = %d, max_dif_tar = %9.8e\n\n", relax_2D.n_x, relax_2D.n_y, max_dif_tar);

    /*Initialize the function values of all grid points*/
    hx = (n_xy.x[n_xy.n_x] - n_xy.x[0]) / n_xy.n_x;  //grid cell width, used in constant and variable terms of relaxation equation
    hy = (n_xy.y[n_xy.n_y] - n_xy.y[0]) / n_xy.n_y;

    x = (float *) malloc((relax_2D.n_x + 1) * sizeof(float));
    y = (float *) malloc((relax_2D.n_y + 1) * sizeof(float));
    f = (float **)malloc((relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= relax_2D.n_x; k++)
    {
        f[k] = (float *)malloc((relax_2D.n_y + 1) * sizeof(float));
    }
    f_star = (float **)malloc((relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= relax_2D.n_x; k++)
    {
        f_star[k] = (float *)malloc((relax_2D.n_y + 1) * sizeof(float));
    }
    f_temp = (float **)malloc((relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= relax_2D.n_x; k++)
    {
        f_temp[k] = (float *)malloc((relax_2D.n_y + 1) * sizeof(float));
    }
    dif = (float **)malloc((relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= relax_2D.n_x; k++)
    {
        dif[k] = (float *)malloc((relax_2D.n_y + 1) * sizeof(float));
    }

    for(i = 0; i <= relax_2D.n_x; i++)             //initialize boundary condition. others initialized to be 0
    {
        for(j = 0; j <= relax_2D.n_y; j++)
        {
            if(i == 0)
                f[i][j] = 0;
            else if(i == relax_2D.n_x)
                f[i][j] = 0;
            else if(j == 0)
                f[i][j] = 0;
            else if(j == relax_2D.n_y)
                f[i][j] = 0;
            else
                f[i][j] = 0;
        }
    }

    /*Iterate until reach the target maximum difference*/
    iter = 0;          //no iterations before calculation
    do
    {
        max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

        for(i = 1; i <= relax_2D.n_x - 1; i++)    //leave boundary unchanged, cause that does not satisfy the equation
        {
            for(j = 1; j <= relax_2D.n_y - 1; j++)
            {
                f_star[i][j] = coef_minus_x * f[i - 1][j] + coef_plus_x * f[i + 1][j] + coef_minus_y * f[i][j - 1] + coef_plus_y * f[i][j + 1] + (*func_ptr)(i, j, hx, hy);
                dif[i][j] = f_star[i][j] - f[i][j];
            }
        }

        for(i = 1; i <= relax_2D.n_x - 1; i++)    //get the maximum difference of these grid points for the present and previous step
        {
            for(j = 1; j <= relax_2D.n_y - 1; j++)
            {
                if(max_dif < fabs(dif[i][j]))
                    max_dif = fabs(dif[i][j]);
            }
        }

        for(i = 1; i <= relax_2D.n_x - 1; i++)        //renew the f(x) with the new f*(x), swap these two arrays
        {
            for(j = 1; j <= relax_2D.n_y - 1; j++)
            {
                f_temp[i][j] = f_star[i][j];
                f_star[i][j] = f[i][j];
                f[i][j] = f_temp[i][j];                     //intermediate process
                //f[i] = f_star[i];                     //direct renew
            }
        }
        //printf("%f\n", f[2]);

        //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);       //add if want to know the intermediate process
        //fprintf(fp, "       x           f(x)       \n");
        //for(i = 0; i <= one_D_n_cell; i++)
            //fprintf(fp, "%9.8e, %9.8e\n", x[i], f[i]);
        //fprintf(fp, "\n");
        iter++;
    } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference

    printf("%d", iter);
    fprintf(fp, "iterations = %d\n\n", iter);
    //printf("%f", max_dif);

    /*Return and save the iteration results*/
    relax_2D.z = (float **)malloc((relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= relax_2D.n_x; k++)
        relax_2D.z[k] = (float *)malloc((relax_2D.n_y + 1) * sizeof(float));
    for(i = 1; i <= relax_2D.n_x - 1; i++)        //renew the f(x) with the new f*(x), swap these two arrays
    {
        for(j = 1; j <= relax_2D.n_y - 1; j++)
        {
                relax_2D.z[i][j] = f[i][j];
        }
    }

    //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);
    fprintf(fp, "       x             y           f(x)       \n");
    for(i = 0; i <= relax_2D.n_x; i++)
    {
        for(j = 0; j <= relax_2D.n_y; j++)
        {
            fprintf(fp, "%9.8e %9.8e %9.8e\n", relax_2D.x[i], relax_2D.y[j], relax_2D.z[i][j]);  //save the final results in data file
            //fprintf(fp, "%9.8e\n", relax_2D.z[i][j]);  //save the final results in data file
        }
    }

    fclose(fp);

    return relax_2D;

    free(x);
    free(f);
    free(f_star);
    free(f_temp);
    free(dif);
};

/*
Two dimensional Gauss-Seidel relaxation method

Two dimensional relaxation method by replacing immediately the iteration result

f_star(x, y) = coef_plus_x * f(x + hx, y) + coef_minus_x * f_star(x - hx, y) + coef_plus_y * f(x, y + hy) + coef_minus_x * f_star(x, y - hy) + g(x, y, a, b)

since the step of x and y may not necessarily be the same.
*/
struct two_dim_grid_data GS_relaxation_2D(float coef_plus_x, float coef_minus_x, float coef_plus_y, float coef_minus_y, float (*func_ptr)(int x_i, int y_j, float hx, float hy))
{
    //void input_two_D_Dir_BC();  //input boundary conditions with different positions
    //void input_two_dim_grid();  //input the number of two dimensional grid cells
    void input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    float hx;       //x width of grid cells
    float hy;       //y width of grid cells
    int i;         //index i of points (x_i, y_j, f(x_i, y_j))
    int j;         //index j of points (x_i, y_j, f(x_i, y_j))
    float *x;      //position x_i
    float *y;      //position y_j
    float **f;      //function values of the prior step
    float f_star;   //function values of the current step, since renew f[i][j] immediately, we needn't storage for new matrix
    float max_dif; //the maximum difference of function values between the current and the prior step of all grid points

    struct two_dim_grid_data gs_relax_2D;  //output result of relaxation method for the 2D boundary value problem
    FILE *fp;                      //output data file for the whole iteration process and the result
    char filename[64];             //filename string

    int iter;          //number of iterations

    /*Input parameters of relaxation method of 2D boundary value problems*/
    //input_two_D_Dir_BC();  //input boundary conditions with different positions   //how
    //input_one_dim_grid();  //input the number of one dimensional grid cells
    input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    /*Create a data file saving the iteration (process and) result*/
    gs_relax_2D.n_x = n_xy.n_x;
    gs_relax_2D.n_y = n_xy.n_y;
    gs_relax_2D.x = n_xy.x;
    gs_relax_2D.y = n_xy.y;
    strcpy(gs_relax_2D.name, result_name);

    sprintf(filename, "GS_relax_2D_%s.txt", gs_relax_2D.name);
    fp = fopen(filename, "w");
    fprintf(fp, "Iteration result of %s using 2D Gauss-Seidel relaxation method\n\n", gs_relax_2D.name);
    fprintf(fp, "Parameters : \n");
    fprintf(fp, "x \\in [%9.8e, %9.8e], y \\in [%9.8e, %9.8e]\n", gs_relax_2D.x[0], gs_relax_2D.x[gs_relax_2D.n_x], gs_relax_2D.y[0], gs_relax_2D.y[gs_relax_2D.n_y]);
    fprintf(fp, "Boundary value : all 0 \n");  //
    fprintf(fp, "n_x = %d, n_y = %d, max_dif_tar = %9.8e\n\n", gs_relax_2D.n_x, gs_relax_2D.n_y, max_dif_tar);

    /*Initialize the function values of all grid points*/
    hx = (n_xy.x[n_xy.n_x] - n_xy.x[0]) / n_xy.n_x;  //grid cell width, used in constant and variable terms of relaxation equation
    hy = (n_xy.y[n_xy.n_y] - n_xy.y[0]) / n_xy.n_y;

    x = (float *) malloc((gs_relax_2D.n_x + 1) * sizeof(float));
    y = (float *) malloc((gs_relax_2D.n_y + 1) * sizeof(float));
    f = (float **)malloc((gs_relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= gs_relax_2D.n_x; k++)
    {
        f[k] = (float *)malloc((gs_relax_2D.n_y + 1) * sizeof(float));
    }

    for(i = 0; i <= gs_relax_2D.n_x; i++)             //initialize boundary condition. others initialized to be 0
    {
        for(j = 0; j <= gs_relax_2D.n_y; j++)
        {
            if(i == 0)
                f[i][j] = 0;
            else if(i == gs_relax_2D.n_x)
                f[i][j] = 0;
            else if(j == 0)
                f[i][j] = 0;
            else if(j == gs_relax_2D.n_y)
                f[i][j] = 0;
            else
                f[i][j] = 0;
        }
    }

    /*Iterate until reach the target maximum difference*/
    iter = 0;          //no iterations before calculation
    do
    {
        max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

        float omega = 0;
        //float omega = 0.93;
        for(i = 1; i <= gs_relax_2D.n_x - 1; i++)    //leave boundary unchanged, cause that does not satisfy the equation
        {
            for(j = 1; j <= gs_relax_2D.n_y - 1; j++)
            {
                f_star = (coef_minus_x * f[i - 1][j] + coef_plus_x * f[i + 1][j] + coef_minus_y * f[i][j - 1] + coef_plus_y * f[i][j + 1] + (*func_ptr)(i, j, hx, hy)) * (1 + omega) - omega * f[i][j];

                if(max_dif < fabs(f_star - f[i][j]))
                    max_dif = fabs(f_star - f[i][j]);

                f[i][j] = f_star;  //renew immediately
            }
        }

        iter++;
    } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference

    printf("%d", iter);
    fprintf(fp, "iterations = %d\n\n", iter);
    //printf("%f", max_dif);

    /*Return and save the iteration results*/
    gs_relax_2D.z = (float **)malloc((gs_relax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= gs_relax_2D.n_x; k++)
        gs_relax_2D.z[k] = (float *)malloc((gs_relax_2D.n_y + 1) * sizeof(float));
    for(i = 1; i <= gs_relax_2D.n_x - 1; i++)        //renew the f(x) with the new f*(x), swap these two arrays
    {
        for(j = 1; j <= gs_relax_2D.n_y - 1; j++)
        {
                gs_relax_2D.z[i][j] = f[i][j];
        }
    }

    //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);
    fprintf(fp, "       x             y           f(x)       \n");
    for(i = 0; i <= gs_relax_2D.n_x; i++)
    {
        for(j = 0; j <= gs_relax_2D.n_y; j++)
        {
            fprintf(fp, "%9.8e %9.8e %9.8e\n", gs_relax_2D.x[i], gs_relax_2D.y[j], gs_relax_2D.z[i][j]);  //save the final results in data file
            //fprintf(fp, "%9.8e\n", gs_relax_2D.z[i][j]);  //save the final results in data file
        }
    }

    fclose(fp);

    return gs_relax_2D;

    free(x);
    free(y);
    free(f);
};

/*
Two dimensional Gauss-Seidel overrelaxation method

Replace f immediately

f_star(x, y) = (1 + \omega) * (coef_plus_x * f(x + hx, y) + coef_minus_x * f_star(x - hx, y) + coef_plus_y * f(x, y + hy) + coef_minus_x * f_star(x, y - hy) + g(x, y, a, b)) - \omega * f(x, y)

since the step of x and y may not necessarily be the same.

\omega = 0 is the normal relaxation method. Here \omega \in (-1, 1).
For \omega \in (0, 1), it is called overrelaxation, which speeds up iteration.
For \omega \in (-1, 0), it also converges, but won't speed up.

Here we use the golden ratio search method to find the best parameter \omega with the fewest iterations.
*/
struct two_dim_grid_data GS_overrelaxation_2D(float coef_plus_x, float coef_minus_x, float coef_plus_y, float coef_minus_y, float (*func_ptr)(int x_i, int y_j, float hx, float hy))
{
    //void input_two_D_Dir_BC();  //input boundary conditions with different positions
    //void input_two_dim_grid();  //input the number of two dimensional grid cells
    void input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    float hx;       //x width of grid cells
    float hy;       //y width of grid cells
    int i;         //index i of points (x_i, y_j, f(x_i, y_j))
    int j;         //index j of points (x_i, y_j, f(x_i, y_j))
    float *x;      //position x_i
    float *y;      //position y_j
    float **f;      //function values of the prior step
    float f_star;   //function values of the current step, since renew f[i][j] immediately, we needn't storage for new matrix
    float max_dif; //the maximum difference of function values between the current and the prior step of all grid points

    struct two_dim_grid_data gs_overrelax_2D;  //output result of relaxation method for the 2D boundary value problem
    FILE *fp;                      //output data file for the whole iteration process and the result
    char filename[64];             //filename string

    float omega[4];       //overrelaxation parameter for golden search, usually start with [0, 0.999]
    int iter[4];          //number of iterations, initialized to be -1 so that we can know if they have been calculated before
    float delta_omega;     //tolerance for the result of \omega
    int next;             //index for the next golden ratio search parameter

    /*Input parameters of relaxation method of 2D boundary value problems*/
    //input_two_D_Dir_BC();  //input boundary conditions with different positions   //how
    //input_one_dim_grid();  //input the number of one dimensional grid cells
    input_max_dif_iter();  //input the target maximum difference of all grid points between two steps

    /*Input the initial points for golden ratio search for parameter \omega*/
    printf("\nEnter the start points for golden ratio search for the best Gauss-Seidel overrelaxation parameter \\omega : \n");
    printf("omega_lower = ");
    scanf("%f", &omega[0]);
    printf("omega_upper = ");
    scanf("%f", &omega[3]);
    omega[1] = (3 - sqrt(5)) * omega[3] / 2 + (sqrt(5) - 1) * omega[0] / 2;
    omega[2] = (3 - sqrt(5)) * omega[0] / 2 + (sqrt(5) - 1) * omega[3] / 2;
    printf("Enter a tolerance for the result delta_omega = ");
    scanf("%f", &delta_omega);
    printf("\n");

    /*Create a data file saving the iteration (process and) result*/
    gs_overrelax_2D.n_x = n_xy.n_x;
    gs_overrelax_2D.n_y = n_xy.n_y;
    gs_overrelax_2D.x = n_xy.x;
    gs_overrelax_2D.y = n_xy.y;
    strcpy(gs_overrelax_2D.name, result_name);

    sprintf(filename, "GS_overrelax_2D_%s.txt", gs_overrelax_2D.name);
    fp = fopen(filename, "w");
    fprintf(fp, "Iteration result of %s using 2D Gauss-Seidel relaxation method\n\n", gs_overrelax_2D.name);
    fprintf(fp, "Parameters : \n");
    fprintf(fp, "x \\in [%9.8e, %9.8e], y \\in [%9.8e, %9.8e]\n", gs_overrelax_2D.x[0], gs_overrelax_2D.x[gs_overrelax_2D.n_x], gs_overrelax_2D.y[0], gs_overrelax_2D.y[gs_overrelax_2D.n_y]);
    fprintf(fp, "Boundary value : all 0 \n");  //
    fprintf(fp, "n_x = %d, n_y = %d, max_dif_tar = %9.8e\n\n", gs_overrelax_2D.n_x, gs_overrelax_2D.n_y, max_dif_tar);

    /*Initialize the function values of all grid points*/
    hx = (n_xy.x[n_xy.n_x] - n_xy.x[0]) / n_xy.n_x;  //grid cell width, used in constant and variable terms of relaxation equation
    hy = (n_xy.y[n_xy.n_y] - n_xy.y[0]) / n_xy.n_y;

    x = (float *) malloc((gs_overrelax_2D.n_x + 1) * sizeof(float));
    y = (float *) malloc((gs_overrelax_2D.n_y + 1) * sizeof(float));
    f = (float **)malloc((gs_overrelax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= gs_overrelax_2D.n_x; k++)
    {
        f[k] = (float *)malloc((gs_overrelax_2D.n_y + 1) * sizeof(float));
    }

    /*Initialize iterations for golden ratio points*/
    fprintf(fp, "   omega    iterations\n");
    printf("   omega    iterations\n");
    for(next = 0; next <= 3; next++)
    {
        for(i = 0; i <= gs_overrelax_2D.n_x; i++)             //initialize boundary condition. others initialized to be 0
        {
            for(j = 0; j <= gs_overrelax_2D.n_y; j++)
            {
                if(i == 0)
                    f[i][j] = 0;
                else if(i == gs_overrelax_2D.n_x)
                    f[i][j] = 0;
                else if(j == 0)
                    f[i][j] = 0;
                else if(j == gs_overrelax_2D.n_y)
                    f[i][j] = 0;
                else
                    f[i][j] = 0;
            }
        }

        /*Iterate until reach the target maximum difference*/
        iter[next] = 0;          //no iterations before calculation
        do
        {
            max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

            for(i = 1; i <= gs_overrelax_2D.n_x - 1; i++)    //leave boundary unchanged, cause that does not satisfy the equation
            {
                for(j = 1; j <= gs_overrelax_2D.n_y - 1; j++)
                {
                    f_star = (coef_minus_x * f[i - 1][j] + coef_plus_x * f[i + 1][j] + coef_minus_y * f[i][j - 1] + coef_plus_y * f[i][j + 1] + (*func_ptr)(i, j, hx, hy)) * (1 + omega[next]) - omega[next] * f[i][j];

                    if(max_dif < fabs(f_star - f[i][j]))
                        max_dif = fabs(f_star - f[i][j]);

                    f[i][j] = f_star;  //renew immediately
                }
            }

            iter[next] = iter[next] + 1;
        } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference

        fprintf(fp, "%7.6e %d\n", omega[next], iter[next]);
        printf("%7.6e %d\n", omega[next], iter[next]);
    }

    /*Golden ratio sear for the best omega and iteration*/
    for(;omega[3] - omega[0] >= delta_omega;)
    {
        next = 0;
        for(int k = 0; k <= 3; k++)
        {
            if(iter[k] < iter[next])
                next = k;
        }

        if(next <= 1)
        {
            next = 1;
            omega[3] = omega[2];
            omega[2] = omega[1];
            omega[1] = (3 - sqrt(5)) * omega[3] / 2 + (sqrt(5) - 1) * omega[0] / 2;
            iter[3] = iter[2];
            iter[2] = iter[1];
        }
        else
        {
            next = 2;
            omega[0] = omega[1];
            omega[1] = omega[2];
            omega[2] = (3 - sqrt(5)) * omega[0] / 2 + (sqrt(5) - 1) * omega[3] / 2;
            iter[0] = iter[1];
            iter[1] = iter[2];
        }

        for(i = 0; i <= gs_overrelax_2D.n_x; i++)             //initialize boundary condition. others initialized to be 0
        {
            for(j = 0; j <= gs_overrelax_2D.n_y; j++)
            {
                if(i == 0)
                    f[i][j] = 0;
                else if(i == gs_overrelax_2D.n_x)
                    f[i][j] = 0;
                else if(j == 0)
                    f[i][j] = 0;
                else if(j == gs_overrelax_2D.n_y)
                    f[i][j] = 0;
                else
                    f[i][j] = 0;
            }
        }

        /*Iterate until reach the target maximum difference*/
        iter[next] = 0;          //no iterations before calculation
        do
        {
            max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

            for(i = 1; i <= gs_overrelax_2D.n_x - 1; i++)    //leave boundary unchanged, cause that does not satisfy the equation
            {
                for(j = 1; j <= gs_overrelax_2D.n_y - 1; j++)
                {
                    f_star = (coef_minus_x * f[i - 1][j] + coef_plus_x * f[i + 1][j] + coef_minus_y * f[i][j - 1] + coef_plus_y * f[i][j + 1] + (*func_ptr)(i, j, hx, hy)) * (1 + omega[next]) - omega[next] * f[i][j];

                    if(max_dif < fabs(f_star - f[i][j]))
                        max_dif = fabs(f_star - f[i][j]);

                    f[i][j] = f_star;  //renew immediately
                }
            }

            iter[next] = iter[next] + 1;
        } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference

        fprintf(fp, "%7.6e %d\n", omega[next], iter[next]);
        printf("%7.6e %d\n", omega[next], iter[next]);
    }

    /*Find the best parameter omega and calculate again the iteration result of that value*/
    int next_temp = next;
    next = 0;
    for(int k = 0; k <= 3; k++)
    {
        if(iter[k] < iter[next])
            next = k;
    }

    if(next != next_temp)
    {
        for(i = 0; i <= gs_overrelax_2D.n_x; i++)             //initialize boundary condition. others initialized to be 0
        {
            for(j = 0; j <= gs_overrelax_2D.n_y; j++)
            {
                if(i == 0)
                    f[i][j] = 0;
                else if(i == gs_overrelax_2D.n_x)
                    f[i][j] = 0;
                else if(j == 0)
                    f[i][j] = 0;
                else if(j == gs_overrelax_2D.n_y)
                    f[i][j] = 0;
                else
                    f[i][j] = 0;
            }
        }

        /*Iterate until reach the target maximum difference*/
        iter[next] = 0;          //no iterations before calculation
        do
        {
            max_dif = 0;                              //every loop need to reset the maximum for this loop since it converges, or the program won't stop

            for(i = 1; i <= gs_overrelax_2D.n_x - 1; i++)    //leave boundary unchanged, cause that does not satisfy the equation
            {
                for(j = 1; j <= gs_overrelax_2D.n_y - 1; j++)
                {
                    f_star = (coef_minus_x * f[i - 1][j] + coef_plus_x * f[i + 1][j] + coef_minus_y * f[i][j - 1] + coef_plus_y * f[i][j + 1] + (*func_ptr)(i, j, hx, hy)) * (1 + omega[next]) - omega[next] * f[i][j];

                    if(max_dif < fabs(f_star - f[i][j]))
                        max_dif = fabs(f_star - f[i][j]);

                    f[i][j] = f_star;  //renew immediately
                }
            }

            iter[next] = iter[next] + 1;
        } while (max_dif >= max_dif_tar);             //stop when reaching the target maximum difference
    }

    printf("\n\\omega_b = %7.6e, iteration = %d\n", omega[next], iter[next]);
    fprintf(fp, "\n\\omega_b = %7.6e, iteration =  %d\n\n",omega[next], iter[next]);
    //printf("%f", max_dif);

    /*Return and save the iteration results*/
    gs_overrelax_2D.z = (float **)malloc((gs_overrelax_2D.n_x + 1) * sizeof(float *));
    for (int k = 0; k <= gs_overrelax_2D.n_x; k++)
        gs_overrelax_2D.z[k] = (float *)malloc((gs_overrelax_2D.n_y + 1) * sizeof(float));
    for(i = 1; i <= gs_overrelax_2D.n_x - 1; i++)        //renew the f(x) with the new f*(x), swap these two arrays
    {
        for(j = 1; j <= gs_overrelax_2D.n_y - 1; j++)
        {
            gs_overrelax_2D.z[i][j] = f[i][j];
        }
    }

    //fprintf(fp, "max_dif = %9.8e\n\n", max_dif);
    fprintf(fp, "       x             y           f(x)       \n");
    for(i = 0; i <= gs_overrelax_2D.n_x; i++)
    {
        for(j = 0; j <= gs_overrelax_2D.n_y; j++)
        {
            fprintf(fp, "%9.8e %9.8e %9.8e\n", gs_overrelax_2D.x[i], gs_overrelax_2D.y[j], gs_overrelax_2D.z[i][j]);  //save the final results in data file
            //fprintf(fp, "%9.8e\n", gs_overrelax_2D.z[i][j]);  //save the final results in data file
        }
    }

    fclose(fp);

    return gs_overrelax_2D;

    free(x);
    free(y);
    free(f);
};

/*
Input 1D Dirichlet boundary condition

f_left = f(x_left)
f_right = f(x_right)
*/
void input_one_D_Dir_BC()
{
    do
    {
        printf("Enter the left boundary condition : \n");
        printf("x_left = ");
        scanf("%f", &one_D_Dir_BC.x_left);
        printf("f_left = ");
        scanf("%f", &one_D_Dir_BC.f_left);
        printf("Enter the right boundary condition : \n");
        printf("x_right = ");
        scanf("%f", &one_D_Dir_BC.x_right);
        printf("f_right = ");
        scanf("%f", &one_D_Dir_BC.f_right);

        if (one_D_Dir_BC.x_left == one_D_Dir_BC.x_right)    //position of two boundaries should be different
            printf("\nPlease enter two different boundary positions x_left and x_right!\n");
        else
            break;
    } while (one_D_Dir_BC.x_left == one_D_Dir_BC.x_right);
};

/*
Input 2D Dirichlet boundary condition with    /////depends, by file or hand or by four functions, given boundary points, calculate by function
*/
void input_two_D_Dir_BC()
{

}

/*
Input one dimensional grid cell number
*/
void input_one_dim_grid()
{
    float n_cell_temp;  //template of the input number of grid cells to determine whether it is an integer

    /*Input an integer larger than one for the number of grid cells*/
    do
    {
        printf("Enter the number of grid cells : n_cell = ");
        scanf("%f", &n_cell_temp);
        one_D_n_cell = (int) n_cell_temp;

        if((one_D_n_cell - n_cell_temp != 0) || (one_D_n_cell <= 1))  //an integer over 1 is needed for positions apart from boundaries
        {
            printf("\nPlease enter an integer n_cell >= 2 !\n");
            one_D_n_cell = 0;
        }
        else
            break;
    } while (one_D_n_cell = 0);
}

/*
Input n dimensional grid cell numbers
*/
void input_multi_dim_grid()
{
    int dim;            //the dimension of independent variables
    float *n_cell;      //number of cells for each dimension
    float n_cell_temp;  //template of the input number of grid cells to determine whether it is an integer

    /*Input dimension of independent variables*/
    do
    {
        printf("Enter the dimension of independent variables : dim = ");
        scanf("%f", &n_cell_temp);
        dim = (int) n_cell_temp;

        if((dim - n_cell_temp != 0) || (dim <= 0))  //dim should be positive integer
        {
            printf("\nPlease enter an integer n_cell >= 2 !\n");
            dim = 0;
        }
        else
            break;
    } while (dim = 0);

    /*Input an integer larger than one for the number of grid cells*/
    n_cell = (float *)malloc(dim * sizeof(float));
    for(int i = 1; i <= dim; i++)
    {
        do
        {
            printf("Enter the number of grid cells : n_cell_%d = ", i);
            scanf("%f", &n_cell_temp);
            n_cell[i - 1]  = (int) n_cell_temp;

            if((n_cell[i - 1] - n_cell_temp != 0) || (n_cell[i - 1] <= 1))  //an integer over 1 is needed for positions apart from boundaries
            {
                printf("\nPlease enter an integer n_cell_%d >= 2 !\n", i);
                n_cell[i - 1] = 0;
            }
            else
                break;
        } while (n_cell[i - 1] = 0);
    }
}

/*
Input the target positive maximum difference of all grid points between two steps of iteration
*/
void input_max_dif_iter()
{
    do
    {
        printf("Enter the target maximum difference between two adjacent steps : max_dif_tar = ");
        scanf("%f", &max_dif_tar);

        if (max_dif_tar <= 0)       //max_dif should be non-positive
            printf("\nPlease enter a positive max_dif_tar!\n");
        else
            break;
    } while (max_dif_tar <= 0);
}

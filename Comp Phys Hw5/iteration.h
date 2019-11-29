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

#include "multi_dim_data.h"
//#include "data_file.h"

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

f_star(x) = coef_plus * f(x - h) + coef_minus * f(x + h) + g(x, h)

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

    /*Input parameters of relaxation method of 1D boundary value problems*/
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
            f_star[i] = coef_minus * f[i - 1] + coef_plus * f[i + 1] + func_ptr(x[i], h);
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
        fprintf(fp, "%9.8e, %9.8e\n", relax_1D.x[i], relax_1D.y[i]);  //save the final results in data file

    fclose(fp);

    return relax_1D;

    free(x);
    free(f);
    free(f_star);
    free(f_temp);
    free(dif);
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

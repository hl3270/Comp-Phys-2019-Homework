/*
multi_dim_data.h

Structures and functions of multi-dimension data.

header files required:
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/

//variables
/*
Structure of n+1 2D data points
*/
struct two_dim_data
{
    int n;             //the number of steps n of data points x_0, ..., x_n
    float *x;          //position x_j, j = 0, ..., n
    float *y;          //function value y_j, j = 0, ..., n

    char name[64];     //name of these data
};

/*
Structure of n+1 3D data points
*/
struct three_dim_data
{
    int n;             //the number of cells
    float *x;          //position x_j, j = 0, ..., n
    float *y;          //position y_j, j = 0, ..., n
    float *z;          //function value z_j, j = 0, ..., n

    char name[64];     //name of these data
};

/*
Structure of 2D grid data
*/
struct two_dim_grid_data
{
    int n_x;           //the number of x cells
    int n_y;           //the number of y cells
    float *x;          //position x_i, i = 0, ..., n_x
    float *y;          //position y_j, j = 0, ..., n_y
    float **z;         //function value z

    char name[64];     //name of these data
};

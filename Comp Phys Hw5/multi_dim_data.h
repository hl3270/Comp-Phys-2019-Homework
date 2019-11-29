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
    int n;             //the number of data points (x_j, y_j)
    float *x;          //position x_j, j = 0, ..., n
    float *y;          //position y_j, j = 0, ..., n
    float *z;          //function value z_j, j = 0, ..., n

    char name[64];     //name of these data
};

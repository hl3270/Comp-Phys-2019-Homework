/*
particle_density.h

Different methods of calculating particle density in multi-dimensional grids

1) Nearest Grid Point (NGP)
Sum up all particles in a grid to calculate the density of the center of the grid

2) Cloud in cell(CIC)
Draw all the grids, and also draw a grid with its center being the position of a particle, and its shape the same as a cell.
Sum up all the overlapped size of the grid cell and all particle grids ti calculate the density of the center of the grid

Header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <string.h>  //string routines
#include <math.h>    //math routines
#include <malloc.h>  //memory allocation routines
*/

#include "multi_dim_data.h"
#include "data_file.h"

#define DATA_FILE_SUFFIX ".txt"

//variables
struct two_dim_grid_data n_xy;

//functions
//struct one_dim_grid_value nearest_grid_point_one_d(){};

//struct two_dim_grid_value nearest_grid_point_two_d(){};

//struct one_dim_grid_value cloud_in_cell_one_d(){};

/*
Calculate the particle density of two dimension grids using cloud-in-cell method, of a square area

Draw all the grids, and also draw a grid with its center being the position of a particle, and its shape the same as a cell.
Sum up all the overlapped size of the grid cell and all particle grids ti calculate the density of the center of the grid
*/
struct two_dim_grid_data cloud_in_cell_two_d()
{
    float x_lower;   //lower boundary of x
    float x_upper;   //upper boundary of x
    float y_lower;   //lower boundary of y
    float y_upper;   //upper boundary of y
    int n_x;         //the number of x grids
    int n_y;         //the number of y grids

    float delta_x;   //x width of a grid
    float delta_y;   //y width of a grid
    float *x;        //center positions x of grids
    float *y;        //center positions y of grids
    float **value;   //value at each grid point

    float xp;        //position x of a particle
    float yp;        //position y of a particle
    int i;           //index i of the cell (i, j) the particle is in
    int j;           //index j of the cell (i, j) the particle is in
    int n_particle;  //number of particles

    struct two_dim_grid_data twoD_gridv;
    FILE *fp;        //data file to save the density of grid points
    char filename[64];  //data file name

    /*Input parameters of grids*/
    printf("Enter the lower boundary of x : x_lower = ");
    scanf("%f", &x_lower);
    printf("Enter the upper boundary of x : x_upper = ");
    scanf("%f", &x_upper);
    printf("Enter the number of grid cells of x : n_x = ");
    scanf("%d", &n_x);
    printf("Enter the lower boundary of y : y_lower = ");
    scanf("%f", &y_lower);
    printf("Enter the upper boundary of y : y_upper = ");
    scanf("%f", &y_upper);
    printf("Enter the number of grid cells of y : n_y = ");
    scanf("%d", &n_y);
    delta_x = (x_upper - x_lower) / n_x;
    delta_y = (y_upper - y_lower) / n_y;

    //printf("%f", delta_x);

    /*Input data of positions of particles*/
    printf("\nInput particle positions by file\n");
    fn = open_data_file(DATA_FILE_SUFFIX);
    //sprintf(filename, "particles.dat.txt");
    fp = fopen(fn.filename, "r");
    n_particle = row_file(fn.filename);
    //printf("%d", n_particle);

    x = (float *) malloc ((n_x + 1) * sizeof(float));
    y = (float *) malloc ((n_y + 1) * sizeof(float));
    value = (float **)malloc((n_x + 1) * sizeof(float *));  //memory allocation for value[n_x + 1]s[n_y + 1]
    for (int k = 0; k <= n_x; k++)
        value[k] = (float *)malloc((n_y + 1) * sizeof(float));

    for(int k = 0; k <= n_x; k++)
        x[k] = x_upper * k / n_x + x_lower * (n_x - k) / n_x;
    for(int k = 0; k <= n_y; k++)
        y[k] = y_upper * k / n_y + y_lower * (n_y - k) / n_y;
    for(int k = 1; k <= n_particle; k++)
    {
        fscanf(fp, "%f %f\n", &xp, &yp);

        //i = (int) (xp - x_lower) * n_x / (x_upper - x_lower);  //not use delta_x for higher accuracy
        //j = (int) (yp - y_lower) * n_y / (y_upper - y_lower);  //not use delta_y for higher accuracy
        //xp = xp - (x_upper - x_lower) * i  / n_x;              //x position of the particle in the cell with origin (x_i, y)j)
        //yp = yp - (y_upper - y_lower) * j  / n_y;              //x position of the particle in the cell with origin (x_i, y)j)
        i = (int) (xp - x_lower) / delta_x;  //not use delta_x for higher accuracy
        j = (int) (yp - y_lower) / delta_y;  //not use delta_y for higher accuracy
        xp = xp - delta_x * i;              //x position of the particle in the cell with origin (x_i, y)j)
        yp = yp - delta_y * j;              //x position of the particle in the cell with origin (x_i, y)j)

        value[i][j] += (delta_x - xp) * (delta_y - yp) / (delta_x * delta_x * delta_y * delta_y);  //add weighted number to grid points
        value[i + 1][j] += xp * (delta_y - yp) / (delta_x * delta_x * delta_y * delta_y);
        value[i][j + 1] += (delta_x - xp) * yp / (delta_x * delta_x * delta_y * delta_y);
        value[i + 1][j + 1] += xp * yp / (delta_x * delta_x * delta_y * delta_y);
    }
    fclose(fp);

    twoD_gridv.n_x = n_x;
    twoD_gridv.n_y = n_y;
    twoD_gridv.x = x;
    twoD_gridv.y = y;
    twoD_gridv.z = value;
    strcpy(twoD_gridv.name, "particle density");

    sprintf(fn.filename, "particle_density_two_dim_CIC.txt");
    fp = fopen(fn.filename, "w");
    fprintf(fp, "Particle density of %s using cloud-in-cell method\n\n");
    fprintf(fp, "x \\in [%f, %f], y \\in [%f, %f]\n", x_lower, x_upper, y_lower, y_upper);
    fprintf(fp, "n_x = %d, n_y = %d, Delta_x = %f, Delta_y = %f\n\n", n_x, n_y, delta_x, delta_y);
    fprintf(fp, "        x             y            density  \n");

    for(i = 0; i <= n_x; i++)
    {
        for(j = 0; j <= n_y; j++)
        {
            fprintf(fp, "%7.6e %7.6e %7.6e\n", x[i], y[j], value[i][j]);
            //fprintf(fp, "%7.6e\n", value[i][j]);
        }
    }

    fclose(fp);

    return twoD_gridv;

    free(x);
    free(y);
    free(value);
};

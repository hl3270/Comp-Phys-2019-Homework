/*
cubic_spline.h

Calculate the cubic spline extrapolation function y = f(x) from a set of points (x_j, y_j), j = 0, 1, 2, ..., n,
and the boundary condition: (1) f'(x_0) and f'(x_n) or (2) f"(x_0) and f"(x_n).
In each interval [x_j, x_j+1], f(x) = a_j * (x - x_j) ^ 3 + b_j * (x - x_j) ^ 2 + c_j * (x - x_j) + d_j, d_j = f(x_j)

Assume x_j is monotonically increasing, or sort x_j and y_j in the increasing order of x_j before calculating
h_j = x_j+1 - x_j, m_j = f"(x_j)
Assume derivatives are continuous at x_j (j != 0, n), we have n-1 equations with n+1 unknowns m_j
h_j-1 * m_j-1 + 2 * (h_j-1 + h_j) * m_j + h_j * m_j+1 = 6 * (((y_j+1 - y_j) / h_j) - ((y_j - y_j-1) / h_j-1))

To determine f(x) exactly, two additional boundary conditions commonly used are:
(1) given f'(x_0), f'(x_j)
(2) given f"(x_0), f"(x_n), specifically, when f"(x_0) = f"(x_n) = 0, it is called natural boundary conditions
(3) periodic boundary condition
(4) not-a-knot condition : [x_0, x_2] uses the same function, so does [x_n-2, x_n]

Since the coefficient matrix is strictly diagonally dominant, after applying boundary conditions,
we can always use Thomas algorithm (tridiagonal LU decomposition) to calculate m_j = f"(x_j)

Then we can use LU decomposition to calculate every f"(x_j) to get the whole expression of f(x)
f(x | [x_j , x_j+1])
= m_j+1 * (x - x_j) ^ 3 / (6 * h_j) - m_j * (x - x_j+1) ^ 3 / (6 * h_j)
  + (6 * y_j+1 - m_j+1 * h_j ^ 2) * (x - x_j) / (6 * h_j)
  - (6 * y_j - m_j * h_j ^ 2) * (x - x_j+1) / (6 * h_j)
= (m_j+1 - m_j) * (x - x_j) ^ 3 / (6 * h_j) + m_j * (x - x_j) ^ 2 / 2
  + ((y_j+1 - y_j) / h_j - (2 * m_j * h_j + m_j+1 * h_j) / 6) * (x - x_j) + y_j

f'(x | [x_j , x_j+1])
= m_j+1 * (x - x_j) ^ 2 / (2 * h_j) - m_j * (x - x_j+1) ^ 2 / (2 * h_j)
  + (y_j+1 - y_j) / h_j - (m_j+1 - m_j) * h_j / 6

header files required:
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
*/

#include "linear_system.h"

//float *Thomas_LU(struct tridiag_matrix trid_mat, float * d);  //LU decomposition solution of tridiagonal linear systems

//variables
/*
Structure of n+1 2D data points
*/
struct two_dim_data
{
    int n;             //the number of steps n of data points x_0, ..., x_n
    float *x;          //position x_j, j = 0, ..., n
    float *y;          //function value y_j, j = 0, ..., n

    char name[64];        //name of these data
};

/*
Structure of output of cubic spline
*/
struct cspline
{
    int n;             //the number of steps n of data points x_0, ..., x_n
    float *x;          //position x_j, j = 0, ..., n
    float *y;          //function value y_j, j = 0, ..., n

    float *h;          //step size h_j

    float *m;          //second derivatives at each position m_j = f"(x_j)

    float *a;          //coefficient of (x - x_j) ^ 3 term in f(x | [x_j, x_j+1])
    float *b;          //coefficient of (x - x_j) ^ 2 term in f(x | [x_j, x_j+1])
    float *c;          //coefficient of (x - x_j) term in f(x | [x_j, x_j+1])

    char name[64];     //name of data
};

//functions
/*
Cubic spline with boundary conditions of two of f'(x_0), f'(x_n), f"(x_0), f"(x_n)
*/

/*
Cubic spline with boundary conditions of first derivatives f'(x_0) & f'(x_n)
*/

/*
Cubic spline with boundary conditions of second derivatives f"(x_0) = m_0, f"(x_n) = m_n, and different step sizes h_j

Delete the 1st and (n+1)th columns of the coefficient matrix
by combining the boundary condition to the constant terms d_j (j = 1,..., n):
d_1 = d_1 - h_0 * f"(x_0)
d_n-1 = d_n-1 - h_n * f"(x_n)

Calculate, save, and return all the other second derivatives m_j's, j = 1, ..., n-1
*/
struct cspline cspline_scd_bc(struct two_dim_data f, float m_0, float m_n)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store

    int j;                  //index j of x_j, y_j, h_j, m_j, etc.
    struct cspline cs;      //cubic spline output

    struct tridiag_matrix coef;  //tridiagonal coefficient (n-1)*(n-1) matrix
    float d[f.n - 1];            //constant term vector
    float *mx;            //solution vector, i.e. second derivatives m_j = f"(x_j), except the two boundary conditions; use pointer since calculated from function

    /*Copy data from input for later use*/
    cs.n = f.n;
    strcpy(cs.name, f.name);
    cs.x = (float *) malloc((cs.n + 1) * sizeof(float));
    cs.y = (float *) malloc((cs.n + 1) * sizeof(float));
    for(j = 0; j <= cs.n; j++)
    {
        cs.x[j] = f.x[j];
        cs.y[j] = f.y[j];
        //printf("%7.6e\n", cs.y[j]);
    }

    /*initiate coefficient matrix and constant vector*/
    cs.h = (float *) malloc(cs.n * sizeof(float));
    for(j = 0; j <= cs.n - 1; j++)
        cs.h[j] = cs.x[j + 1] - cs.x[j];

    coef.n = f.n - 1;
    coef.a = (float *) malloc((coef.n - 1) * sizeof(float));
    coef.c = (float *) malloc((coef.n - 1) * sizeof(float));
    coef.b = (float *) malloc(coef.n * sizeof(float));
    for(j = 0; j <= coef.n - 2; j++)
    {
        coef.a[j] = cs.h[j + 1];
        coef.c[j] = cs.h[j + 1];
        coef.b[j] = 2 * (cs.h[j] + cs.h[j + 1]);
    }
    coef.b[coef.n - 1] = 2 * (cs.h[coef.n - 1] + cs.h[coef.n]);
    //printf("%7.6e\n", coef.a[coef.n - 2]);

    d[0] = 6 * ((cs.y[2] - cs.y[1]) / cs.h[1] - (cs.y[1] - cs.y[0]) / cs.h[0]) - cs.h[0] * m_0;  //combined with boundary condition
    d[cs.n - 2] = 6 * ((cs.y[cs.n] - cs.y[cs.n - 1]) / cs.h[cs.n - 1] - (cs.y[cs.n - 1] - cs.y[cs.n - 2]) / cs.h[cs.n - 2]) - cs.h[f.n - 1] * m_n;
    for(j = 1; j <= cs.n - 3; j++)
        d[j] = 6 * ((cs.y[j + 2] - cs.y[j + 1]) / cs.h[j + 1] - (cs.y[j + 1] - cs.y[j]) / cs.h[j]);
        //printf("%7.6e\n", d[cs.n - 4]);

    /*Calculate m_1, ..., m_n-1*/
    mx = Thomas_LU(coef, d);     //mx[j] = f"(x_j+1) = m_j+1

    /*Calculate and store the coefficients a_j b_j c_j ,in f(x|[x_j, x_j+1]) = a_j*(x-x_j)^3 + b_j*(x-x_j)^2 + c_j*(x-x_j) + y_j*/
    cs.m = (float *) malloc((cs.n + 1) * sizeof(float));
    cs.m[0] = m_0;
    cs.m[cs.n] = m_n;
    for(j = 1; j <= cs.n - 1; j++)
        cs.m[j] = mx[j - 1];

    cs.a = (float *) malloc((cs.n) * sizeof(float));
    cs.b = (float *) malloc((cs.n) * sizeof(float));
    cs.c = (float *) malloc((cs.n) * sizeof(float));
    for(j = 0; j <= cs.n - 1; j++)
    {
        cs.a[j] = (cs.m[j + 1] - cs.m[j]) / (6 * cs.h[j]);
        cs.b[j] = cs.m[j] / 2;
        cs.c[j] = (cs.y[j + 1] - cs.y[j]) / cs.h[j] - (2 * cs.m[j] + cs.m[j + 1]) * cs.h[j] / 6;
    }

    /*Save and return the results*/
    sprintf(filename, "cubic_spline_of_%s.txt", cs.name);
    fp = fopen(filename, "w");
    fprintf(fp, "Cubic spline of %s with boundary conditions : f\"(x_0) = %7.6e, f\"(x_n) = %7.6e\n\n", cs.name, m_0, m_n);
    fprintf(fp, "f(x | [x_j , x_j+1]) \n= m_j+1 * (x - x_j) ^ 3 / (6 * h_j) - m_j * (x - x_j+1) ^ 3 / (6 * h_j) \n");
    fprintf(fp, "  + (6 * y_j+1 - m_j+1 * h_j ^ 2) * (x - x_j) / (6 * h_j) \n  - (6 * y_j - m_j * h_j ^ 2) * (x - x_j+1) / (6 * h_j)\n");
    fprintf(fp, "= (m_j+1 - m_j) * (x - x_j) ^ 3 / (6 * h_j) + m_j * (x - x_j) ^ 2 / 2 \n");
    fprintf(fp, "  + ((y_j+1 - y_j) / h_j - (2 * m_j * h_j + m_j+1 * h_j) / 6) * (x - x_j) + y_j\n\n");
    fprintf(fp, "    j       x_j            y_j      h_j=x_j+1-x_j   m_j=f\"(x_j)        a_j            b_j            c_j      \n");

    for(j = 0; j <= cs.n - 1; j++)
        fprintf(fp, "%5d %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n", j, cs.x[j], cs.y[j], cs.h[j], cs.m[j], cs.a[j], cs.b[j], cs.c[j]);
    fprintf(fp, "%5d %7.6e %7.6e                %7.6e", j, cs.x[cs.n], cs.y[cs.n], cs.m[cs.n]);

    fclose(fp);

    return cs;
}

/*
Cubic spline with natural boundary conditions f"(x_0) = f"(x_n) = 0

Directly delete the 1st and (n+1)th columns of the coefficient matrix

Calculate, save, and return all the other second derivatives m_j's
*/
struct cspline cspline_ntr_bc(struct two_dim_data f)
{
    return cspline_scd_bc(f, 0, 0);
}

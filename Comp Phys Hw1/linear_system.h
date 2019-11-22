/*
linear_system.h

Solving linear systems

Required header files:
#include <stdio.h>   //standard input/output library
#include <stdlib.h>  //standard library
#include <math.h>    //standard math library
#include <malloc.h>  //memory allocation routines
*/

/*
Structure of a tridiagonal (n * n) matrix
b_0 c_0
a_0 b_1 c_1
    a_1 b_2 c_2
        .   .   .
         .   .   .
          .   .   .
           .   .   .
            a_n-3 b_n-2 c_n-2
                  a_n-2 b_n-1
stored in sub-diagonal, diagonal, and super-diagonal arrays a[n-1], b[n], c[n-1]
*/
struct tridiag_matrix
{
    int n;           //the index n of the (n * n) matrix
    float *a;        //sub-diagonal array a[n - 1]
    float *b;        //diagonal array b[n]
    float *c;        //super-diagonal array c[n - 1]
};

/*
Gauss elimination
*/

/*
LU decomposition
*/

/*
Chasing method

LU decomposition for general banded matrix
*/

/*
Thomas algorithm

LU decomposition for tridiagonal matrices
a[n - 1], b[n], c[n - 1] are sub-diagonal, diagonal, and super-diagonal arrays of tridiagonal matrix,
d[n] is the constant term array of the non-homogeneous linear system

The following algorithm is based on matrix A = LU_0, i.e. U matrix is a unit upper triangular matrix
NB: the sub-diagonal terms of L matrix are the same as that of A matrix : alpha_i = a_i

Process: Assume U_0x = y;
         Calculate l[n] diagonal terms of L matrix, u[n - 1] super-diagonal terms of U_0 matrix, and vector y[n]
         Calculate x[n] with back substitution

For tridiagonal matrices, the condition is : det(A) != 0      ///////add pivoting if not diagonally dominant
Or: |b[0]| > |c[0]| > 0;
    |b[k]| > |a[k]| + |c[k]|, a[k] * c[k] != 0 (1 <= k <= n - 2);
    |b[n - 1]| > |a[n - 2]| > 0
*/
float *Thomas_LU(struct tridiag_matrix trid_mat, float *d)    //add : if n != size of d, dimension unequal return false
{
    int i;
    float l;                   //diagonal array of L matrix   //needn't be array if don't need the expression of
    float u[trid_mat.n - 1];   //super-diagonal array of U_0 matrix
    float y[trid_mat.n];       //medium vector in solving LU decomposition linear system
    float *x;                  //solution vector

    /*Calculate L matrix (won't store), U_0 matrix, and y vector*/
    l = trid_mat.b[0];         //l_1 = b_1, here index j is row index
    y[0] = d[0] / l;           //y_1 = d_1 / l_1
    for(i = 0; i <= trid_mat.n - 2; i++)
    {
        u[i] = trid_mat.c[i] / l;                             //u_j = c_j / l_j
        l = trid_mat.b[i + 1] - trid_mat.a[i] * u[i];         //l_j+1 = b_j+1 - a_j+1 * u_j
        y[i + 1] = (d[i + 1] - trid_mat.a[i] * y[i]) / l;
    }

    /*Calculate solution vector x by back substitution*/
    x = (float *) malloc(trid_mat.n * sizeof(float));
    x[trid_mat.n - 1] = y[trid_mat.n - 1];
    for(i = trid_mat.n - 2; i >= 0; i--)
        x[i] = y[i] - u[i] * x[i + 1];

    return x;
}

//matrix pivoting

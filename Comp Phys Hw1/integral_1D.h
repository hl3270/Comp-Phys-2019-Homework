/*
integral_1D.h

1D definite integrals over [a,b]

Only return the integration result

Below are header files required
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <math.h>    //math routines
*/

float a;                       //lower bound of the integral
float b;                       //upper bound of the integral
int N_0;                       //number of bins
float eps_tol;                 //user defined tolerance for the relative error of the integral
int max_loop;                  //user defined maximum number for the loops to be operated

/*
Integration with midpoint rule
*/
float midpoint(float a, float b, float (*func_ptr)(float), int N_0, float eps_tol, int max_loop)
{
    int n;                     //number of bins
    int i;                     //loop number
    float temp;                //template for the weighted sum over function values at sample points
    int j;                     //sample position index
    float x;                   //sample position x_j
    float integral[max_loop];  //integrals for loop i and i+1
    float e[max_loop - 1];     //absolute error of integrals for i and i+1
    float eps[max_loop - 1];   //relative error of integrals for i and i+1
    int i_max = max_loop;      //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;                //epsilon index template used to find the smallest epsilon

    n = N_0;                   //to separate from the N widely used in global definition

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                             //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= n; j++)                              //weighted sum over all function values at sample points
        {
            x = a + (j + 1 / 2) * (b - a) / n;                //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
            integral[0] = temp * (b - a) / n;                 //integrate using midpoint rule
        else
        {
            integral[i - 1] = temp * (b - a) / n;             //integrate using midpoint rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
        }

        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        n = n * 2;             //double the number of bins
    }

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                      //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                              //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Integration with trapezoid rule
*/
float trapezoid(float a, float b, float (*func_ptr)(float), int N_0, float eps_tol, int max_loop)
{
    int n;                     //number of bins
    int i;                     //loop number
    float temp;                //template for the weighted sum over function values at sample points
    int j;                     //sample position index
    float x;                   //sample position x_j
    float integral[max_loop];  //integrals for loop i and i+1
    float e[max_loop - 1];     //absolute error of integrals for i and i+1
    float eps[max_loop - 1];   //relative error of integrals for i and i+1
    int i_max = max_loop;      //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;                //epsilon index template used to find the smallest epsilon

    n = N_0;                   //to separate from the N widely used in global definition

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                             //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= n; j++)                              //weighted sum over all function values at sample points
        {
            x = a + j * (b - a) / n;                          //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
            integral[0] = temp * (b - a) / n - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * n);            //integrate using midpoint rule
        else
        {
            integral[i - 1] = temp * (b - a) / n - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * n);             //integrate using midpoint rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
        }

        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        n = n * 2;             //double the number of bins
    }

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                      //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                              //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Integration with Simpson's rule
*/
float Simpson(float a, float b, float (*func_ptr)(float), int N_0, float eps_tol, int max_loop)
{
    int n;                     //number of bins
    int i;                     //loop number
    float temp;                //template for the weighted sum over function values at sample points
    int j;                     //sample position index
    float x;                   //sample position x_j
    float integral[max_loop];  //integrals for loop i and i+1
    float e[max_loop - 1];     //absolute error of integrals for i and i+1
    float eps[max_loop - 1];   //relative error of integrals for i and i+1
    int i_max = max_loop;      //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;                //epsilon index template used to find the smallest epsilon

    n = N_0;                   //to separate from the N widely used in global definition

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                             //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= n; j++)                              //weighted sum over all function values at sample points
        {
            x = a + j * (b - a) / n;                          //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            if (j % 2 == 1)
                temp = temp + 4 * (*func_ptr)(x);             //for odd index j, the weight of the function value in the sum is 4
            else
                temp = temp + 2 * (*func_ptr)(x);             //for even index j, the weight of the function value in the sum is 2
        }

        if(i == 1)
            integral[0] = temp * (b - a) / (3 * n) - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (3 * n);            //integrate using midpoint rule
        else
        {
            integral[i - 1] = temp * (b - a) / (3 * n) - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (3 * n);             //integrate using midpoint rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
        }

        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        n = n * 2;             //double the number of bins
    }

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                      //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                              //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Integration with midpoint method
*/
float Romberg(float a, float b, float (*func_ptr)(float), int N_0, float eps_tol, int max_loop)
{
    int n;                     //number of bins
    int i;                     //index i for R_i,j , i.e. loop number for trapezoid rule R_i,1 = I_i
    int j;                     //index j for R_i,j
    int k;                     //sample position index of trapezoid rule
    float x;                   //sample position x_k of trapezoid rule
    float temp;                //template for the weighted sum over function values at sample points
    float R[max_loop * (max_loop + 1) / 2];  //store R_i,j in an array in the order of R_1,1 R_2,1 R_2,2 R_3,1 R_3,2 R_3,3 etc.
    float e[max_loop - 1];     //absolute error of integrals R_i,i
    float eps[max_loop - 1];   //relative error of integrals R_i,i
    int i_max = max_loop;      //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;                //epsilon index template used to find the smallest epsilon

    n = N_0;                   //to separate from the N widely used in global definition

    for (i = 1; i <= max_loop; i++)
    {
        for (j = 1; j <= i; j++)
        {
            if (j == 1)        //R_i,1 = I_i of trapezoid rule
            {
                temp = 0;
                for (k = 0; k <= n; k++)
                {
                    x = a + k * (b - a) / n;
                    temp = temp + (*func_ptr)(x);
                }
                R[(i * (i - 1) / 2)] = temp * (b - a) / n - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * n);
            }
            else               //calculate R_i,j (j > 1)
                R[(i * (i - 1) / 2) + j - 1] = R[(i * (i - 1) / 2) + j - 2] +(R[(i * (i - 1) / 2) + j - 2] - R[((i - 1) * (i - 2) / 2) + j - 2]) / ((4 ^ (j - 1) - 1));

            if (j > 1)         //break condition 1: stop calculating when R_i,j-1 - R_i-1,j-1 reaches machine precision to avoid overflow, and return the maximum of loops really calculated
            {
                if (fabs(R[(i * (i - 1) / 2) + j - 2] - R[((i - 1) * (i - 2) / 2) + j - 2]) <= (R[((i - 1) * (i - 2) / 2) + j - 2] * ((4 ^ (j - 1) - 1)) * (2 ^ (- 23))))    //if the term (R_i,j-1 - R_i,j-2) / (4^j-1 - 1) reaches machine precision, then stop //still need fixing, won't stop immediately when out put overflow
                {
                    i_max = i;
                    break;
                }
            }

            e[i - 2] = R[(i * (i + 1) / 2) - 1] - R[(i * (i - 1) / 2) - 1];
            eps[i - 2] = e[i - 2] / R[(i * (i - 1) / 2) - 1];
        }

        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))  //break condition 2: when reaching the user defined convergence condition
        {
            i_max = i;
            break;
        }

        n = n * 2;
    }

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                      //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return R[(((i_temp + 2) * (i_temp + 3)) / 2) - 1];        //eps[i] is related to R_i+2,i+2 i.e. R[((i_temp + 2) * (i_temp + 3) / 2) - 1]
}

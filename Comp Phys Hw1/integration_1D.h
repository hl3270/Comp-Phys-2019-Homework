/*
integration_1D.h

Including methods of calculating definite integrals over [a,b]

Data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon

Below are header files required
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <math.h>    //math routines
#include <stdbool.h> //standard boolean library
#include <string.h>  //string library
*/


//variables
float a;                    //lower bound of the integral
float b;                    //upper bound of the integral
float (*func_ptr)(float x); //function pointer
int N;                      //number of bins
float eps_tol;              //user defined tolerance for the relative error of the integral
int max_loop;               //user defined maximum number for the loops to be operated

int int_method;        //choose which method to use for integration
float N_temp;            //template of N to check if have input a positive integer
float max_loop_temp;     //template of max_loop to check if have input a positive integer


//functions
void midpoint_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral_1, integral_2;   //integrals for loop i and i+1
    float e_1, e_2;                 //absolute error of integrals for i and i+1
    float eps_1 = eps_tol, eps_2;   //relative error of integrals for i and i+1

    sprintf(filename, "mid_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using midpoint rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N - 1; j++)                     //weighted sum over all function values at sample points
        {
            x = - a - (j + 1 / 2) * (b - a) / N;         //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
        {
            integral_1 = temp * (b - a) / N;             //integrate using midpoint rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral_1);
        }
        else
        {
            integral_2 = temp * (b - a) / N;             //integrate using midpoint rule

            if(i == 2)
            {
                e_1 = integral_2 - integral_1;           //initiate absolute error of the integral
                eps_1 = e_1 / integral_1;                //initiate relative error of the integral
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_1, eps_1);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = integral_2 - integral_1;            //absolute error of the integral for loop i + 1
                eps_2 = e_2 / integral_1;                 //relative error of the integral for loop i + 1
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_2, eps_2);
                integral_1 = integral_2;                  //for the next loop i+1, integral_2 = I_i+1 in this loop i should take the place of integral_1 = I_i, since the new i = i+1
                e_1 = e_2;                                //similar to above, e_2 for this loop should take the place of e_1 in the next loop
                eps_1 = eps_2;                            //similar to above, eps_2 for this loop should take the place of eps_1 in the next loop
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        N = N * 2;                      //double the number of bins
    }

    fclose(fp);
}

void trapezoid_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral_1, integral_2;   //integrals for loop i and i+1
    float e_1, e_2;         //absolute error of integrals for i and i+1
    float eps_1 = eps_tol, eps_2;     //relative error of integrals for i and i+1

    sprintf(filename, "tpz_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using trapezoid rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N - 1; j++)                     //weighted sum over all function values at sample points
        {
            x = - a - j * (b - a) / N;                   //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
        {
            integral_1 = temp * (b - a) / N + ((*func_ptr)(- a) + (*func_ptr)(- b)) * (b - a) / (2 * N);  //integrate using trapezoid rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral_1);
        }
        else
        {
            integral_2 = temp * (b - a) / N + ((*func_ptr)(- a) + (*func_ptr)(- b)) * (b - a) / (2 * N);  //integrate using trapezoid rule

            if(i == 2)
            {
                e_1 = integral_2 - integral_1;           //initiate absolute error of the integral
                eps_1 = e_1 / integral_1;                //initiate relative error of the integral
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_1, eps_1);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = integral_2 - integral_1;            //absolute error of the integral for loop i + 1
                eps_2 = e_2 / integral_1;                 //relative error of the integral for loop i + 1
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_2, eps_2);
                integral_1 = integral_2;                  //for the next loop i+1, integral_2 = I_i+1 in this loop i should take the place of integral_1 = I_i, since the new i = i+1
                e_1 = e_2;                                //similar to above, e_2 for this loop should take the place of e_1 in the next loop
                eps_1 = eps_2;                            //similar to above, eps_2 for this loop should take the place of eps_1 in the next loop
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        N = N * 2;                      //double the number of bins
    }

    fclose(fp);
}

void Simpson_s_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral_1, integral_2;   //integrals for loop i and i+1
    float e_1, e_2;         //absolute error of integrals for i and i+1
    float eps_1 = eps_tol, eps_2;     //relative error of integrals for i and i+1

    sprintf(filename, "Sps_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using Simpson's rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N - 1; j++)                     //weighted sum over all function values at sample points
        {
            x = - a - j * (b - a) / N;                   //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            if (j % 2 == 1)
                temp = temp + 4 * (*func_ptr)(x);           //for odd index j, the weight of the function value in the sum is 4
            else
                temp = temp + 2 * (*func_ptr)(x);           //for even index j, the weight of the function value in the sum is 2
        }

        if(i == 1)
        {
            integral_1 = temp * (b - a) / (3 * N) + ((*func_ptr)(- a) + (*func_ptr)(- b)) * (b - a) / (3 * N);  //integrate using Simpson's rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral_1);
        }
        else
        {
            integral_2 = temp * (b - a) / (3 * N) + ((*func_ptr)(- a) + (*func_ptr)(- b)) * (b - a) / (3 * N);  //integrate using Simpson's rule

            if(i == 2)
            {
                e_1 = integral_2 - integral_1;           //initiate absolute error of the integral
                eps_1 = e_1 / integral_1;                //initiate relative error of the integral
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_1, eps_1);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = integral_2 - integral_1;            //absolute error of the integral for loop i + 1
                eps_2 = e_2 / integral_1;                 //relative error of the integral for loop i + 1
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_2, eps_2);
                integral_1 = integral_2;                  //for the next loop i+1, integral_2 = I_i+1 in this loop i should take the place of integral_1 = I_i, since the new i = i+1
                e_1 = e_2;                                //similar to above, e_2 for this loop should take the place of e_1 in the next loop
                eps_1 = eps_2;                            //similar to above, eps_2 for this loop should take the place of eps_1 in the next loop
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        N = N * 2;                      //double the number of bins
    }

    fclose(fp);
}

//unfinished romberg
void Romberg_integration_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral_1, integral_2;   //integrals for loop i and i+1
    float e_1, e_2;                 //absolute error of integrals for i and i+1
    float eps_1 = eps_tol, eps_2;   //relative error of integrals for i and i+1

    sprintf(filename, "mid_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using midpoint rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N - 1; j++)                     //weighted sum over all function values at sample points
        {
            x = - a - (j + 1 / 2) * (b - a) / N;         //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
        {
            integral_1 = temp * (b - a) / N;             //integrate using midpoint rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral_1);
        }
        else
        {
            integral_2 = temp * (b - a) / N;             //integrate using midpoint rule

            if(i == 2)
            {
                e_1 = integral_2 - integral_1;           //initiate absolute error of the integral
                eps_1 = e_1 / integral_1;                //initiate relative error of the integral
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_1, eps_1);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = integral_2 - integral_1;            //absolute error of the integral for loop i + 1
                eps_2 = e_2 / integral_1;                 //relative error of the integral for loop i + 1
                fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral_2, e_2, eps_2);
                integral_1 = integral_2;                  //for the next loop i+1, integral_2 = I_i+1 in this loop i should take the place of integral_1 = I_i, since the new i = i+1
                e_1 = e_2;                                //similar to above, e_2 for this loop should take the place of e_1 in the next loop
                eps_1 = eps_2;                            //similar to above, eps_2 for this loop should take the place of eps_1 in the next loop
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        N = N * 2;                      //double the number of bins
    }

    fclose(fp);
}

int choose_int_method()  //choose which method to use for integration
{
    int i_method;      //integration method index. 1: midpoint rule; 2: trapezoid rule; 3: Simpson's rule; 4: Romberg integration

    do
    {
        printf("Which method do you want to choose for integration? \n");
        printf("For midpoint rule, type 1; \nFor trapezoid rule, type 2; \nFor Simpson's rule, type 3; \nFor Romberg Integration, type 4. \n");  //choose algorithm for differentiation
        scanf("%d",&i_method);
        printf("\n");

        switch (i_method)
        {
        case 1:
            return i_method;
        case 2:
            return i_method;
        case 3:
            return i_method;
        default:
            printf("Wrong method choice!\n");
            i_method = 0;
        }
    } while (i_method == 0);
}

float input_int_lower_bound()     //input lower bound a of the integral
{
    float a;

    printf("Enter the lower bound of the integral  a = ");
    scanf("%f", &a);
    printf("\n");

    return a;
}

float input_int_upper_bound()     //input upper bound b of the integral
{
    float b;

    printf("Enter the upper bound of the integral  b = ");
    scanf("%f", &b);
    printf("\n");

    return b;
}

int input_bins()       //input the initial number of bins N
{
    int bins;
    float bins_temp;

    do
    {
        printf("Enter the initial number of bins N = ");
        scanf("%f", &bins_temp);
        printf("\n");
        bins = (int) bins_temp;

        if ((bins + 1 - bins_temp != 1) || (bins_temp <= 0))        //find if the input N is a positive integer
        {
            printf("Enter a positive integer!\n");
            bins = 0;
        }
        else
        {
            switch (int_method)
            {
            case 1:
                return bins;
            case 2:
                return bins;
            case 3:
                {
                    if (bins % 2 != 0)                         //find if the input N for Simpson's rule is a positive even integer
                    {
                        printf("Enter an positive even integer for Simpson's rule!\n");
                        bins = 0;
                    }
                    else
                        return bins;
                }
            default:
                break;
            }
        }
    } while (bins == 0);
}

float input_epsilon_tolerance()   //input tolerance for relative error
{
    float eps_tol;

    do
    {
        printf("Enter a tolerance for the relative error eps_tol = ");
        scanf("%f", &eps_tol);
        printf("\n");

        switch(eps_tol == 0)
        {
        case true:
            printf("Enter a non-zero number!\n\n");
        default:
            return eps_tol;
        }
    } while (eps_tol == 0);
}

int input_max_of_loop()  //input a user defined maximum for the number of loops of calculation
{
    int max_loop;
    float max_loop_temp;

    do
    {
        printf("Enter in case a maximum number for loops of calculation max_loop = ");
        scanf("%f", &max_loop_temp);
        printf("\n");
        max_loop = (int) max_loop_temp;

        switch (max_loop + 1 - max_loop_temp != 1 || max_loop_temp <= 0)
        {
        case true:
            printf("Enter a positive integer!\n");
        default:
            return max_loop;
        }
    } while (max_loop + 1 - max_loop_temp != 1 || max_loop_temp <= 0);
}

void calculate_integral(int int_method, float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)  //calculate integrals in the chosen method with parameters entered
{
    switch (int_method)
    {
    case 1:
        midpoint_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        break;
    case 2:
        trapezoid_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        break;
    case 3:
        Simpson_s_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        break;
    case 4:
        Romberg_integration_1D(a, b, func_ptr, N, eps_tol, max_loop);
        break;
    default:
        break;
    }
}

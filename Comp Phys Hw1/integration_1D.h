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


//variables     ////////insert in one structure?????
float a;                    //lower bound of the integral
float b;                    //upper bound of the integral
float (*func_ptr)(float x); //function pointer
int N;                      //number of bins
float eps_tol;              //user defined tolerance for the relative error of the integral
int max_loop;               //user defined maximum number for the loops to be operated

int int_method;             //choose which method to use for integration


//functions
/*
Calculate integral of f(x) within [a,b] using midpoint rule
*/
float midpoint_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral[max_loop];//integrals for loop i and i+1
    float e[max_loop - 1];  //absolute error of integrals for i and i+1
    float eps[max_loop - 1];//relative error of integrals for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "mid_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using midpoint rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N; j++)                         //weighted sum over all function values at sample points
        {
            x = a + (j + 1 / 2) * (b - a) / N;           //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
        {
            integral[0] = temp * (b - a) / N;            //integrate using midpoint rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral[0]);
        }
        else
        {
            integral[i - 1] = temp * (b - a) / N;             //integrate using midpoint rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
            fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        N = N * 2;           //double the number of bins
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                 //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                         //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Calculate integral of f(x) within [a,b] using trapezoid rule
*/
float trapezoid_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral[max_loop];//integrals for loop i and i+1
    float e[max_loop - 1];  //absolute error of integrals for i and i+1
    float eps[max_loop - 1];//relative error of integrals for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "tpz_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using trapezoid rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N; j++)                         //weighted sum over all function values at sample points
        {
            x = a + j * (b - a) / N;                     //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            temp = temp + (*func_ptr)(x);
        }

        if(i == 1)
        {
            integral[0] = temp * (b - a) / N - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * N);  //integrate using trapezoid rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral[0]);
        }
        else
        {
            integral[i - 1] = temp * (b - a) / N - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * N);  //integrate using trapezoid rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
            fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        N = N * 2;             //double the number of bins
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                 //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                         //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Calculate integral of f(x) within [a,b] using Simpson's rule
*/
float Simpson_s_rule_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //loop number
    float temp;             //template for the weighted sum over function values at sample points
    int j;                  //sample position index
    float x;                //sample position x_j
    float integral[max_loop];//integrals for loop i and i+1
    float e[max_loop - 1];  //absolute error of integrals for i and i+1
    float eps[max_loop - 1];//relative error of integrals for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "Sps_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using Simpson's rule\n\n", function_1D, a, b);
    fprintf(fp, " i        N  integral       e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = 0;                                        //initiate the template for sum as 0, before adding any function value
        for (j = 0; j <= N; j++)                         //weighted sum over all function values at sample points
        {
            x = a + j * (b - a) / N;                     //sample position, calculated repeatedly to decrease error of t when N grows large and j -> N
            if (j % 2 == 1)
                temp = temp + 4 * (*func_ptr)(x);        //for odd index j, the weight of the function value in the sum is 4
            else
                temp = temp + 2 * (*func_ptr)(x);        //for even index j, the weight of the function value in the sum is 2
        }

        if(i == 1)
        {
            integral[0] = temp * (b - a) / (3 * N) - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (3 * N);  //integrate using Simpson's rule
            fprintf(fp, "%2d  %7d  %7.6e\n", i, N, integral[0]);
        }
        else
        {
            integral[i - 1] = temp * (b - a) / (3 * N) - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (3 * N);  //integrate using Simpson's rule
            e[i - 2] = integral[i - 1] - integral[i - 2];     //absolute error of the integral
            eps[i - 2] = e[i - 2] / integral[i - 2];          //relative error of the integral
            fprintf(fp, "%2d  %7d  %7.6e  %7.6e  %7.6e\n", i, N, integral[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, bins N, integral, estimated absolute error e, relative error \epsilon
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))
        {
            i_max = i;
            break;
        }

        N = N * 2;             //double the number of bins
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= max_loop - 2; l++)              //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return integral[i_temp + 1];                         //eps[i] is related to I_i+2 i.e. integral[i+1]
}

/*
Calculate integral of f(x) within [a,b] using Romberg integration
*/
float Romberg_integration_1D(float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store integrals, numbers of bins, and estimated errors
    int i;                  //index i for R_i,j , i.e. loop number for trapezoid rule R_i,1 = I_i
    int j;                  //index j for R_i,j
    int k;                  //sample position index of trapezoid rule
    float x;                //sample position x_k of trapezoid rule
    float temp;             //template for the weighted sum over function values at sample points
    float R[max_loop * (max_loop + 1) / 2]; //store R_i,j in an array in the order of R_1,1 R_2,1 R_2,2 R_3,1 R_3,2 R_3,3 etc.
    float e[max_loop - 1];  //absolute error of integrals R_i,i
    float eps[max_loop - 1];//relative error of integrals R_i,i
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "Rbg_I_of_%s_BTW_%7.6e_%7.6e.txt", function_1D, a, b);
    fp = fopen(filename, "w");
    fprintf(fp, "Integrals of %s within [%7.6e,%7.6e] using Romberg integration\n\n", function_1D, a, b);
    fprintf(fp, "i,j        N  integral R_i,j  e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        for (j = 1; j <= i; j++)
        {
            if (j == 1)     //R_i,1 = I_i of trapezoid rule
            {
                temp = 0;
                for (k = 0; k <= N; k++)
                {
                    x = a + k * (b - a) / N;
                    temp = temp + (*func_ptr)(x);
                }
                R[(i * (i - 1) / 2)] = temp * (b - a) / N - ((*func_ptr)(a) + (*func_ptr)(b)) * (b - a) / (2 * N);
            }
            else            //calculate R_i,j (j > 1)
                R[(i * (i - 1) / 2) + j - 1] = R[(i * (i - 1) / 2) + j - 2] +(R[(i * (i - 1) / 2) + j - 2] - R[((i - 1) * (i - 2) / 2) + j - 2]) / ((4 ^ (j - 1) - 1));

            if (j > 1)      //break condition 1: stop calculating when R_i,j-1 - R_i-1,j-1 reaches machine precision to avoid overflow, and return the maximum of loops really calculated
            {
                if (fabs(R[(i * (i - 1) / 2) + j - 2] - R[((i - 1) * (i - 2) / 2) + j - 2]) <= (R[((i - 1) * (i - 2) / 2) + j - 2] * ((4 ^ (j - 1) - 1)) * (2 ^ (- 23))))    //if the term (R_i,j-1 - R_i,j-2) / (4^j-1 - 1) reaches machine precision, then stop //still need fixing, won't stop immediately when out put overflow
                {
                    i_max = i;
                    break;
                }
            }

            if ((j != i) || (i == 1))
                fprintf(fp, "%d,%d  %7d  %7.6e \n", i, j, N, R[(i * (i - 1) / 2) + j - 1]);
            else
            {
                e[i - 2] = R[(i * (i + 1) / 2) - 1] - R[(i * (i - 1) / 2) - 1];
                eps[i - 2] = e[i - 2] / R[(i * (i - 1) / 2) - 1];
                fprintf(fp, "%d,%d  %7d  %7.6e  %7.6e  %7.6e \n", i, j, N, R[(i * (i + 1) / 2) - 1], e[i - 2], eps[i - 2]);
            }
        }

        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))  //break condition 2: when reaching the user defined convergence condition
        {
            i_max = i;
            break;
        }

        N = N * 2;
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                 //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return R[(((i_temp + 2) * (i_temp + 3)) / 2) - 1];   //eps[i] is related to R_i+2,i+2 i.e. R[((i_temp + 2) * (i_temp + 3) / 2) - 1]
}

/*
Choose which method to use for integration
*/
int choose_int_method()
{
    int i_method;                 //integration method index. 1: midpoint rule; 2: trapezoid rule; 3: Simpson's rule; 4: Romberg integration

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
        case 4:
            return i_method;
        default:
            printf("Wrong method choice!\n");
            i_method = 0;
        }
    } while (i_method == 0);
}

/*
Input lower bound a of the integral
*/
float input_int_lower_bound()
{
    float a;

    printf("Enter the lower bound of the integral  a = ");
    scanf("%f", &a);
    printf("\n");

    return a;
}

/*
Input upper bound b of the integral
*/
float input_int_upper_bound()
{
    float b;

    printf("Enter the upper bound of the integral  b = ");
    scanf("%f", &b);
    printf("\n");

    return b;
}

/*
Input the initial number of bins N
*/
int input_bins()
{
    int bins;
    float bins_temp;              //template of N to check if have input a positive integer

    do
    {
        printf("Enter the initial number of bins N = ");
        scanf("%f", &bins_temp);
        printf("\n");
        bins = (int) bins_temp;

        if ((bins + 1 - bins_temp != 1) || (bins_temp <= 0))   //find if the input N is a positive integer
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

/*
Input tolerance for relative error
*/
float input_epsilon_tolerance()
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

/*
Input a user defined maximum for the number of loops of calculation
*/
int input_max_of_loop()
{
    int max_loop;
    float max_loop_temp;          //template of max_loop to check if have input a positive integer

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

/*
Integrate f(x) within [a,b] using different methods with initial number of bins N,
until reaching the user defined relative error tolerance eps_tol or until finishing max_loop loops of calculation,
give the best result and save all the results in files
*/
void calculate_integral(int int_method, float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop)  //calculate integrals in the chosen method with parameters entered
{
    float I;

    switch (int_method)
    {
    case 1:
        I = midpoint_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        printf("I = %7.6e\n", I);
        break;
    case 2:
        I = trapezoid_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        printf("I = %7.6e\n", I);
        break;
    case 3:
        I = Simpson_s_rule_1D(a, b, func_ptr, N, eps_tol, max_loop);
        printf("I = %7.6e\n", I);
        break;
    case 4:
        I = Romberg_integration_1D(a, b, func_ptr, N, eps_tol, max_loop);
        printf("I = %7.6e\n", I);
        break;
    default:
        break;
    }
}

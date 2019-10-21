/*
differentiation_1D.h

Including methods of calculating the derivative of f_1D(x) at a position x

Data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon

Below are header files required:
#include <stdlib.h>  //standard library
#include <stdio.h>   //standard input/output library
#include <math.h>    //math routines
#include <stdbool.h> //standard boolean library
#include <string.h>  //string library
*/


//variables     ////////insert in one structure?????
float x;                    //position x
float h;                    //step size
float (*func_ptr)(float x); //function pointer
float eps_tol;              //user defined tolerance for the relative error of the integral
int max_loop;               //user defined maximum number for the loops to be operated

int dif_method;             //choose differentiation method


//functions
/*
Differentiate using forward-difference algorithm
*/
float forward_difference_1D(float x, float h, float (*func_ptr)(float x), float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives and step sizes
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative[max_loop];//derivatives of the function for h_i and h_i+1
    float e[max_loop - 1];  //absolute error of derivatives for i and i+1
    float eps[max_loop - 1]; //relative error of derivatives for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "fwd_d_of_%s_at_%f.txt", function_1D, x);///////////try input manually filename
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using forward difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative[0] = (((*func_ptr)(x + h) - (*func_ptr)(x)) / h);      //differentiate using forward difference method
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative[0]);
        }
        else
        {
            derivative[i - 1] = (((*func_ptr)(x + h) - (*func_ptr)(x)) / h);  //differentiate using forward difference method
            e[i - 2] = derivative[i - 1] - derivative[i - 2];                 //the absolute error of the derivative
            eps[i - 2] = e[i - 2] / derivative[i - 2];                        //the relative error of the derivative
            fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
        }
        ////////////add: when f(x+h) - f(x) or something reaches machine precision then break
        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));         //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))         //break condition: when reaching the user defined convergence condition
        {
            i_max = i;
            break;
        }

        h = h / 2;          //half the step size
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                                      //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return derivative[i_temp + 1];                                            //eps[i] is related to f_prime_i+2 i.e. derivative[i+1]
}

/*
Differentiate using central-difference algorithm
*/
float central_difference_1D(float x, float h, float (*func_ptr)(float x), float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives and step sizes
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative[max_loop];//derivatives of the function for h_i and h_i+1
    float e[max_loop - 1];  //absolute error of derivatives for i and i+1
    float eps[max_loop - 1]; //relative error of derivatives for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "ctr_d_of_%s_at_%f.txt", function_1D, x);
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using central difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative[0] = (((*func_ptr)(x + h) - (*func_ptr)(x - h)) / (2 * h));      //differentiate using central difference method
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative[0]);
        }
        else
        {
            derivative[i - 1] = (((*func_ptr)(x + h) - (*func_ptr)(x - h)) / (2 * h));  //differentiate using central difference method
            e[i - 2] = derivative[i - 1] - derivative[i - 2];                           //the absolute error of the derivative
            eps[i - 2] = e[i - 2] / derivative[i - 2];                                  //the relative error of the derivative
            fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));                   //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))                   //break condition: when reaching the user defined convergence condition
        {
            i_max = i;
            break;
        }

        h = h / 2;          //half the step size
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                                                //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return derivative[i_temp + 1];                                                      //eps[i] is related to f_prime_i+2 i.e. derivative[i+1]
}

/*
Differentiate using extrapolated-difference algorithm
*/
float extrapolated_difference_1D(float x, float h, float (*func_ptr)(float), float eps_tol, int max_loop)
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives, step sizes, and estimated errors
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative[max_loop];//derivatives of the function for h_i and h_i+1
    float e[max_loop - 1];  //absolute error of derivatives for i and i+1
    float eps[max_loop - 1];//relative error of derivatives for i and i+1
    int i_max = max_loop;   //maximum of loop number really calculated, initialized to be max_loop, which is the maximum case; will change if stop before reaching this value
    int i_temp;             //epsilon index template used to find the smallest epsilon

    sprintf(filename, "etp_d_of_%s_at_%f.txt", function_1D, x);
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using extrapolated difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i <= max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative[0] = ((-(*func_ptr)(x + 2 * h) + 8 * (*func_ptr)(x + h) - 8 * (*func_ptr)(x - h) + (*func_ptr)(x - 2 * h)) / (12 * h));      //differentiate using central difference method
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative[0]);
        }
        else
        {
            derivative[i - 1] = ((-(*func_ptr)(x + 2 * h) + 8 * (*func_ptr)(x + h) - 8 * (*func_ptr)(x - h) + (*func_ptr)(x - 2 * h)) / (12 * h));  //differentiate using central difference method
            e[i - 2] = derivative[i - 1] - derivative[i - 2];                            //the absolute error of the derivative
            eps[i - 2] = e[i - 2] / derivative[i - 2];                                   //the relative error of the derivative
            fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative[i - 1], e[i - 2], eps[i - 2]);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));                    //test
        if ((log(fabs(eps[i - 2])) < log(fabs(eps_tol))) && (i != 1))                    //break condition: when reaching the user defined convergence condition
        {
            i_max = i;
            break;
        }

        h = h / 2;          //half the step size
    }

    fclose(fp);

    i_temp = 0;
    for (int l = 1; l <= i_max - 2; l++)                                                 //find index i when relative error is the smallest
    {
        if (fabs(eps[i_temp]) > fabs(eps[l]))
            i_temp = l;
    }
    return derivative[i_temp + 1];                                                       //eps[i] is related to f_prime_i+2 i.e. derivative[i+1]
}

/*
Choose to use which method above to calculate derivatives
*/
int choose_dif_method()
{
    int d_method;           //differentiate method index. 1: forward-difference; 2: central-difference; 3: extrapolated-difference

    do
    {
        printf("Which differentiation method do you want to choose foe differentiation? \nFor forward-difference algorithm, type 1; \nFor central-difference algorithm, type 2; \nFor extrapolated-difference algorithm, type 3. \n");  //choose algorithm for differentiation
        scanf("%d",&d_method);
        printf("\n");

        switch (d_method)
        {
        case 1:
            return d_method;
        case 2:
            return d_method;
        case 3:
            return d_method;
        default:
            printf("Wrong method choice!\n");
            d_method = 0;
        }
    } while (d_method == 0);
}

/*
Input the position x where you want to calculate    /////// add identification of false input
*/
float input_position()
{
    float x;                //position x

    printf("Enter the position x = ");
    scanf("%f", &x);
    printf("\n");

    return x;
}

/*
Input the initial step size h_0
*/
float input_step_size()
{
    float h;                //step size h

    do
    {
        printf("Set initial step size h = ");
        scanf("%f", &h);
        printf("\n");

        switch(h <= 0)
        {
        case true:
            printf("Enter a positive number!\n");
        default:
            return h;
        }
    } while (h <= 0);
}

/*
Input user defined relative error tolerance to determine convergence
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
    float max_loop_temp;

    do
    {
        printf("Enter in case a maximum for the number of loops of calculation max_loop = ");
        scanf("%f", &max_loop_temp);
        printf("\n");
        max_loop = (int) max_loop_temp;

        switch (max_loop + 1 - max_loop_temp != 1 || max_loop_temp <= 0)  //identify if have input a positive integer
        {
        case true:
            printf("Enter a positive integer!\n");
        default:
            return max_loop;
        }
    } while (max_loop + 1 - max_loop_temp != 1 || max_loop_temp <= 0);
}

/*
Calculate derivative of f(x) using different methods at position x with initial step size h,
until reaching a user defined relative error tolerance eps_tol or until finishing max_loop loops of calculation,
give the best result and store all the results in files
*/
void calculate_derivative_1D(int dif_method, float x, float h, float (*func_ptr)(float), float eps_tol, int max_loop)  //calculate derivatives in the chosen method with parameters entered
{
    float f_prime;

    switch (dif_method)                   //use different methods to calculate derivatives
    {
    case 1:
        f_prime = forward_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        printf("f_prime = %7.6e\n", f_prime);
        break;
    case 2:
        f_prime = central_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        printf("f_prime = %7.6e\n", f_prime);
        break;
    case 3:
        f_prime = extrapolated_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        printf("f_prime = %7.6e\n", f_prime);
        break;
    default:
        break;
    }
}

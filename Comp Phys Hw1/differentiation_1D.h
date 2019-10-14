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


//variables     ////////put in one structure?????
float x;                    //position x
float h;                    //step size
float (*func_ptr)(float x); //function pointer
float eps_tol;              //user defined tolerance for the relative error of the integral
int max_loop;               //user defined maximum number for the loops to be operated

int dif_method;             //choose differentiation method
float max_loop_temp;        //template of max_loop to check if have input a positive integer


//functions
void forward_difference_1D(float x, float h, float (*func_ptr)(float x), float eps_tol, int max_loop)  //differentiate cos(x) using forward-difference algorithm
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives and step sizes
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative_1, derivative_2;  //derivatives of the function for h_i and h_i+1
    float e_1, e_2;                    //absolute error of derivatives for i and i+1
    float eps_1 = eps_tol, eps_2;      //relative error of derivatives for i and i+1

    sprintf(filename, "fwd_d_of_%s_at_%f.txt", function_1D, x);////////////////////
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using forward difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i < max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative_1 = (((*func_ptr)(x + h) - (*func_ptr)(x)) / h);  //initiate the derivative for i
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative_1);
        }
        else
        {
            derivative_2 = (((*func_ptr)(x + h) - (*func_ptr)(x)) / h);  //the derivative for i+1

            if(i == 2)
            {
                e_1 = derivative_2 - derivative_1;                       //initiate the absolute error of the derivative
                eps_1 = e_1 / derivative_1;                              //initiate the relative error of the derivative
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_1, eps_1);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = derivative_2 - derivative_1;
                eps_2 = e_2 / derivative_1;
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_2, eps_2);
                e_1 = e_2;                                               //set e_i for the next i=i+1 as e_i+1 for this loop i
                eps_1 = eps_2;                                           //set eps_i for the next i=i+1 as eps_i+1 for this loop i
                derivative_1 = derivative_2;
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        h = h / 2;          //half the step size
    }

    fclose(fp);
}

void central_difference_1D(float x, float h, float (*func_ptr)(float x), float eps_tol, int max_loop)  //differentiate cos(x) using forward-difference algorithm
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives and step sizes
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative_1, derivative_2;  //derivatives of the function for h_i and h_i+1
    float e_1, e_2;                    //absolute error of derivatives for i and i+1
    float eps_1 = eps_tol, eps_2;      //relative error of derivatives for i and i+1

    sprintf(filename, "ctr_d_of_%s_at_%f.txt", function_1D, x);
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using central difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i < max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative_1 = (((*func_ptr)(x + h) - (*func_ptr)(x - h)) / (2 * h));  //initiate the derivative for i
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative_1);
        }
        else
        {
            derivative_2 = (((*func_ptr)(x + h) - (*func_ptr)(x - h)) / (2 * h));  //the derivative for i+1

            if(i == 2)
            {
                e_1 = derivative_2 - derivative_1;                                 //initiate the absolute error of the derivative
                eps_1 = e_1 / derivative_1;                                        //initiate the relative error of the derivative
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_1, eps_1);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = derivative_2 - derivative_1;
                eps_2 = e_2 / derivative_1;
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_2, eps_2);
                e_1 = e_2;                                                         //set e_i for the next i=i+1 as e_i+1 for this loop i
                eps_1 = eps_2;                                                     //set eps_i for the next i=i+1 as eps_i+1 for this loop i
                derivative_1 = derivative_2;
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        h = h / 2;          //half the step size
    }

    fclose(fp);
}

void extrapolated_difference_1D(float x, float h, float (*func_ptr)(float), float eps_tol, int max_loop)  //differentiate cos(x) using forward-difference algorithm
{
    char filename[64];      //filename string
    FILE *fp;               //pointer to file to store derivatives, step sizes, and estimated errors
    int i;                  //loop number
    float temp;             //template used in decreasing the error between h and the difference between x & x+h
    float derivative_1, derivative_2;  //derivatives of the function for h_i and h_i+1
    float e_1, e_2;                    //absolute error of derivatives for i and i+1
    float eps_1 = eps_tol, eps_2;      //relative error of derivatives for i and i+1

    sprintf(filename, "etp_d_of_%s_at_%f.txt", function_1D, x);
    fp = fopen(filename, "w");
    fprintf(fp, "Derivatives of %s at %f using extrapolated difference method\n\n", function_1D, x);
    fprintf(fp, " i  h              derivative      e               epsilon \n");

    for (i = 1; i < max_loop; i++)
    {
        temp = x + h;       //to decrease error of the difference between x & x+h
        h = temp - x;

        if(i == 1)
        {
            derivative_1 = (((*func_ptr)(x + 2 * h) + 8 * (*func_ptr)(x + h) - 8 * (*func_ptr)(x - h) - (*func_ptr)(x - 2 * h)) / (12 * h));  //initiate the derivative for i
            fprintf(fp, "%2d  %7.6e  %7.6e\n", i, h, derivative_1);
        }
        else
        {
            derivative_2 = (((*func_ptr)(x + 2 * h) + 8 * (*func_ptr)(x + h) - 8 * (*func_ptr)(x - h) - (*func_ptr)(x - 2 * h)) / (12 * h));  //the derivative for i+1

            if(i == 2)
            {
                e_1 = derivative_2 - derivative_1;                                 //initiate the absolute error of the derivative
                eps_1 = e_1 / derivative_1;                                        //initiate the relative error of the derivative
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_1, eps_1);  //data are written in the form: loop i, step size h, derivative, estimated absolute error e, relative error \epsilon
            }
            else
            {
                e_2 = derivative_2 - derivative_1;
                eps_2 = e_2 / derivative_1;
                fprintf(fp, "%2d  %7.6e  %7.6e  %7.6e  %7.6e\n", i, h, derivative_2, e_2, eps_2);
                e_1 = e_2;                                                         //set e_i for the next i=i+1 as e_i+1 for this loop i
                eps_1 = eps_2;                                                     //set eps_i for the next i=i+1 as eps_i+1 for this loop i
                derivative_1 = derivative_2;
            }
        }

        //printf("%lf, %lf\n", log(fabs(eps_1)), log(fabs(eps_tol)));  //test
        if (log(fabs(eps_1)) < log(fabs(eps_tol)))
            break;

        h = h / 2;          //half the step size
    }

    fclose(fp);
}

int choose_dif_method()     //choose to use which method above to calculate derivatives
{
    int d_method;             //differentiate method index. 1: forward-difference; 2: central-difference; 3: extrapolated-difference

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

float input_position()      //input the position x where you want to calculate    /////// add identification of false input
{
    float x;                //position x

    printf("Enter the position x = ");
    scanf("%f", &x);
    printf("\n");

    return x;
}

float input_step_size()     //input the initial step size h_0
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

float input_epsilon_tolerance()  //input user defined relative error tolerance to determine convergence
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

int input_max_of_loop()     //input a user defined maximum for the number of loops of calculation
{
    int max_loop;
    float max_loop_temp;

    do
    {
        printf("Enter in case a maximum for the number of loops of calculation max_loop = ");
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

void calculate_derivative_1D(int dif_method, float x, float h, float (*func_ptr)(float), float eps_tol, int max_loop)  //calculate derivatives in the chosen method with parameters entered
{
    switch (dif_method)                   //use different methods to calculate derivatives
    {
    case 1:
        forward_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        break;
    case 2:
        central_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        break;
    case 3:
        extrapolated_difference_1D(x, h, func_ptr, eps_tol, max_loop);
        break;
    default:
        break;
    }
}

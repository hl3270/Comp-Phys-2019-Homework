/* Hw_1-2.c */

/*
Integrate exp(-t) over range t = [0,1] for single precision
using midpoint rule, trapezoid rule, and Simpson's rule.

Compare the relative error \epsilon for these three methods
with log-log plot of \epsilon as a function of numbers of bins N,
and make sure to make N large enough to see the effects of roundoff error.

And explain what you see in the plot.
*/

#include <stdlib.h>       //standard library
#include <stdio.h>        //input/output library
#include <math.h>         //math routines
#include <stdbool.h>      //string boolean library
#include <string.h>       //string library

#include "function_input.h"  //function input header file
#include "integration_1D.h"  //integration header file

void main()
{
    int exit;                //choose whether to exit after a calculation

    //functions
    int choose_function();   //choose which function to calculate
    int choose_int_method(); //choose which method to use for integration
    float f_1D(float x);     //functions to calculate
    char *f_1D_expression(int function_choice);  //expression of the 1D function, used in filename and inside the file
    float input_int_lower_bound();    //input lower bound a of the integral
    float input_int_upper_bound();    //input upper bound b of the integral
    int input_bins();        //input the initial number of bins N_0
    float input_epsilon_tolerance();  //input tolerance for relative error
    int input_max_of_loop(); //input a user defined maximum for the number of loops of calculation
    void calculate_integral(int int_method, float a, float b, float (*func_ptr)(float), int N, float eps_tol, int max_loop);  //calculate integrals in the chosen method with parameters entered

    printf("Integrate exp(-t) over [a,b] \n\n");

    beginning:

        function_choice = choose_function();  //choose which function to calculate
        int_method = choose_int_method();     //choose a method of differentiation

        do
        {
            a = input_int_lower_bound();      //input lower bound a of the integral
            b = input_int_upper_bound();      //input upper bound b of the integral

            if (a == b)                       //input again if a = b
                printf("Please enter different numbers for lower and upper bounds!\n");
        } while (a == b);

        N = input_bins();                     //input the initial number of bins N
        eps_tol = input_epsilon_tolerance();  //input tolerance for relative error
        max_loop = input_max_of_loop();       //input a maximum for loops
        func_ptr = f_1D;                      //choose function
        function_1D = f_1D_expression(function_choice);  //to put the expression of functions into the filename and the file
        calculate_integral(int_method, a, b, func_ptr, N, eps_tol, max_loop);

        printf("\nTo continue, type 1; \nTo exit, press any other key.\n");      //choose to continue or exit
        scanf("%d", &exit);
        if (exit==1)
        {
            printf("\n");
            goto beginning;
        }
        else exit;
}

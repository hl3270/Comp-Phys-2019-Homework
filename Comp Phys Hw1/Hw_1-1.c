/* Hw_1-1.c */

/*
Differentiate the functions cos(x) and exp(x) at x = 0.1, 10,
using single precision forward-, central-, and extrapolated difference algorithms.

Check whether the scaling of relative error \epsilon,
and the number of significant digits obtained agrees with simple estimates,
with log-log plot of \epsilon vs step size h.

And clarify the regimes in the plot where truncation error and roundoff error manifest themselves.
*/

#include <stdlib.h>          //standard library
#include <stdio.h>           //input/output library
#include <math.h>            //math routines
#include <stdbool.h>         //standard boolean library
#include <string.h>          //string library

#include "function_input.h"  //function input header file
#include "differentiation_1D.h" //differentiation header file

void main()
{
    int exit;                //choose whether to exit after a calculation

    //functions
    //char *function_input();
    int choose_function();
    int choose_dif_method(); //choose which method to use for differentiation
    float f_1D(float x);     //functions to calculate
    char *f_1D_expression(int function_choice);  //expression of the 1D function, used in filename and inside the file
    float input_position();  //input position x
    float input_step_size(); //input step size h
    float input_epsilon_tolerance();  //input tolerance for relative error
    int input_max_of_loop(); //input a user defined maximum for the number of loops of calculation
    void calculate_derivative_1D(int dif_method, float x, float h, float (*func_ptr)(float), float eps_tol, int max_loop);  //calculate derivatives with parameters entered

    printf("Calculate derivatives of cos(x) or exp(x) \n\n");

    beginning:
        //function_1D = function_input();
        function_choice = choose_function();  //choose which function to calculate
        dif_method = choose_dif_method();     //choose a method of differentiation
        x = input_position();                 //input position x
        h = input_step_size();                //input step size h
        eps_tol = input_epsilon_tolerance();  //input tolerance for relative error
        max_loop = input_max_of_loop();       //input a maximum for loops
        func_ptr = f_1D;                      //choose function
        function_1D = f_1D_expression(function_choice);  //to put the expression of functions into the filename and the file
        calculate_derivative_1D(dif_method, x, h, func_ptr, eps_tol, max_loop);  //calculate derivatives with parameters entered

        printf("\nTo continue, type 1; \nTo exit, press any other key.\n");  //choose to continue or exit
        scanf("%d", &exit);
        if (exit==1)
        {
            printf("\n");
            goto beginning;
        }
        else exit;
}

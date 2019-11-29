/*Hw_5-1.c*/

/*
Exercise 9.7 (Computational Physics, Mark Newman)
The relaxation method for ordinary differential equations

There is no reason why the relaxation method must be restricted to the solution of differential equations with two or more independent variables.
It can also be applied to those with one independent variable, i.e., to ordinary differential equations.
In this context, as with partial differential equations, it is a technique for solving boundary value problems,
which are less common with ordinary differential equations but do occur - we discussed them in Section 8.6.

Consider the problem we looked at in Example 8.8 on page 390, in which a ball of mass m = 1 kg is thrown
from height x = 0 into the air and lands back at x = 0 ten seconds later.
The problem is to calculate the trajectory of the ball, but we cannot do it using initial value methods
like the ordinary Runge-Kutta method because we are not told the initial velocity of the ball.
One approach to finding a solution is the shooting method of Section 8.6.1. Another is the relaxation method.

Ignoring friction effects, the trajectory is the solution of the ordinary differential equation

d^2 x / d t^2 = - g

where g is the acceleration due to gravity.

a) Replacing the second derivative in this equation with its finite-difference approximation, Eq.(5.109),

f"(x) \approx (f(x + h) - 2 * f(x) + f(x - h)) / h ^ 2

derive a relaxation-method equation for solving this problem on a time-like "grid" of points with separation h.

b) Taking the boundary conditions to be that x = 0 at t = 0 and t = 10, write a program to solve
for the height of the ball as a function of time using the relaxation method with 100 points
and make a plot of the result from t = 0 to t = 10. Run the relaxation method until the answers change
by 10^(-6) or less at every point on each step.

Note that, unlike the shooting method, the relaxation method does not give us the initial value of the velocity needed
to achieve the required solution. It gives us only the solution itself, although one could get an approximation
to the initial velocity by calculating a numerical derivative of the solution at time t = 0.
On balance, however, the relaxation method for ordinary differential equations is most useful when one wants to know
the details of the solution itself, but not the initial conditions needed to achieve it.
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines

#include "iteration.h"  //iteration methods header file

#define GE (9.80665)    //the gravity of earth
#define RESULT_NAME "trajectory"  //name of the iteration result

void main()
{
    struct two_dim_data relaxation_1D(float coef_plus, float coef_minus, float (*func_ptr)(float x, float h));  //relaxation method of one dimensional boundary value problem
    float func_x_h(float x, float h);   //constant term in the relaxation equation. for more functions, can add a variable to choose function

    struct two_dim_data trajectory;     //iteration result

    beginning:
        strcpy(result_name, RESULT_NAME);
        printf("Calculate trajectory due to gravity\n\n");
        trajectory = relaxation_1D(pow(2, -1), pow(2, -1), func_x_h);  //input conditions, parameters, and calculate the relaxation result
        printf("\n");

        for(int i = 0; i <= one_D_n_cell; i++)
            printf("%9.8e %9.8e\n", trajectory.x[i], trajectory.y[i]); //one dimensional function result of relaxation method

        free(trajectory.x);
        free(trajectory.y);

        printf("\nTo continue, type 1; \nTo exit, press any other key.\n");  //choose to continue or exit
        scanf("%d", &exit);
        if (exit==1)
        {
            printf("\n");
            goto beginning;
        }
        else exit;
}

/*
Constant term in relaxation equation

f_star(x) = coef_plus * f(x - h) + coef_minus * f(x + h) + g(x, h)
*/
float func_x_h(float x, float h)
{
    return (GE * h * h / 2);
}

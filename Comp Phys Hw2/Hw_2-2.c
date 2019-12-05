/*Hw_2-2.c*/

/*
Exercise 6.13 (Computational Physics, Mark Newman)
Wien's displacement constant

Planck's radiation law tells us that the intensity of radiation per unit area and per unit wavelength \lambda
from a black body at temperature T is

I(\lambda) = 2 * \pi * h * c ^ 2 * \lambda ^ (-5) / (exp (h * c / (\lambda * kB * T)) - 1)

where h is Planck's constant, c is the speed of light, and k_B is Boltzmann's constant.

a) Show by differentiating that the wavelength \lambda at which the emitted radiation is strongest is the solution of the equation

5 * exp(- h * c / (\lambda * k_B * T) + h * c / (\lambda * k_B * T) - 5 = 0

Make the substitution x = h * c / (\lambda * k_B * T) and hence show that the wavelength of maximum radiation obeys
the Wien displacement law : \lambda = b / T, where the so-called Wien displacement constant is
b = h * c / (k_B * x), and x is the solution to the nonlinear equation 5 * exp(-x) + x - 5 = 0.

b) Write a program to solve this equation to an accuracy of \epsilon = 10^(-6) using the binary search method,
and hence find a value for the displacement constant.

c) The displacement law is the basis for the method of optical pyrometry, a method for measuring the temperatures
of objects by observing the color of the thermal radiation they emit.
The method is commonly used to estimate  the surface temperatures of astronomical bodies, such as the Sun.
The wavelength peak in the Sun's emitted radiation falls at \lambda = 502 nm.
From the equations above and your value of the displacement constant, estimate the surface temperature of the Sun.
*/

#include <stdlib.h>        //standard library
#include <stdio.h>         //standard input/output library
#include <string.h>        //string routines
#include <math.h>          //math routines

#include "root_finding.h"  //root finding header file

void main()
{
    float f(float x);                //the equation f(x) = 0 pending root finding
    float bin_search(float x_1, float x_2, float (*func_ptr)(float x), float accuracy);  //standard relaxation method

    float x_1;                       //left starting point
    float x_2;                       //right starting point
    float x;                         //independent variable
    float (*func_ptr)(float x);      //function pointer
    float accuracy;                  //target accuracy

    func_ptr = f;
    strcpy(equation_name, "5exp(-x)+x-5");
    printf("Solution to %s using binary search\n\n", equation_name);
    //do
    //{
        //printf("Enter the left starting point x_1 = ");      //enter manually the starting point and the target accuracy
        //scanf("%f", &x_1);
        //printf("Enter the right starting point x_2 = ");
        //scanf("%f", &x_2);
        //if((x_1 == x_2) || ((*func_ptr)(x_1) * (*func_ptr)(x_1) > 0))
            //printf("\nPlease enter another two numbers!\n\n");
    //} while ((x_1 == x_2) || ((*func_ptr)(x_1) * (*func_ptr)(x_1) > 0));
    //printf("\nEnter the target solution accuracy \\epsilon = ");
    //scanf("%f", &accuracy);
    //printf("\n");

    x_1 = 4;                         //predefined left starting point
    x_2 = 6;                         //predefined right starting point
    accuracy = 1e-6;                 //predefined target accuracy
    printf("x_1 = %7.6e, x_2 = %7.6e, accuracy = %7.6e\n\n", x_1, x_2, accuracy);
    printf("      x        %s      \n", equation_name);

    x = bin_search(x_1, x_2, func_ptr, accuracy);
    printf("\nx = %7.6e\n\n", x);
}

float f(float x)
{
    return (5 * exp(- x) + x - 5);
}

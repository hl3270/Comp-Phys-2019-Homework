/*Hw_2-1.c*/

/*
Exercise 6.11 (Computational Physics, Mark Newman)
Overrelaxation

If you did not already do Exercise 6.10, you should do it before this one.
The ordinary relaxation method involves iterating the equation x′ = f(x), starting from an initial guess, until it converges.
As we have seen, this is often a fast and easy way to find solutions to nonlinear equations.
However, it is possible in some cases to make the method work even faster using the technique of overrelaxation.
Suppose our initial guess at the solution of a particular equation is, say, x = 1, and the final, true solution is x = 5.
After the first step of the iterative process, we might then see a value of, say, x = 3.
In the overrelaxation method, we observe this value and note that x is increasing,
then we deliberately overshoot the calculated value, in the hope that this will get us closer to the final solution
- in this case we might pass over x = 3 and go straight to a value of x = 4 perhaps,
which is closer to the final solution of x = 5 and hence should get us to that solution quicker.
The overrelaxation method provides a formula for performing this kind of overshooting in a controlled fashion and often,
though not always, it does get us to our solution faster. In detail, it works as follows.

We can rewrite the equation x′ = f(x) in the form x′ = x + \Delta x, where

\Delta x = x′ - x = f(x) - x.

The overrelaxation method involves iteration of the modified equation

x′ = x + (1 + \omega) \Delta x,

(keeping the definition of \Delta x the same). If the parameter \omega is zero, then this is the same as
the ordinary relaxation method, but for \omega > 0 the method takes the amount \Delta x by which the value of x
would have been changed and changes by a little more. Using \Delta x = f(x) - x, we can also write x′ as

x′ = x + (1 + \omega) * [f (x) - x] = (1 + \omega) f(x) - \omega x,

which is the form in which it is usually written.

For the method to work the value of \omega must be chosen correctly, although there is some wiggle room
- there is an optimal value, but other values close to it will typically also give good results.
Unfortunately, there is no general theory that tells us what the optimal value is.
Usually it is found by trial and error.

a) Derive an equivalent of

\epsilon′ = \epsilon * f′(x*), x* = x′ + \epsilon′

for the overrelaxation method and hence show that the error on x′, the equivalent of

\epsilon′ = (x - x′) / (1 - 1 / f′(x*) \approx (x - x′) / (1 - 1 / f′(x)

, is given by

\epsilon′ \approx (x - x′ / (1 - 1 / ((1 + \omega) * f′(x) - \omega))

b) Consider again the equation x = 1 - e ^ (- c * x) that we solved in Exercise 6.10.
Take the program you wrote for part (a) of that exercise, which solved the equation for the case c = 2,
and modify it to print out the number of iterations it takes to converge to a solution accurate to 10^(-6).

c) Now write a new program (or modify the previous one) to solve the same equation x = 1 - e ^ (- c * x) for c = 2,
again to an accuracy of 10^(-6), but this time using overrelaxation.
Have your program print out the answers it finds along with the number of iterations it took to find them.
Experiment with different values of \omega to see how fast you can get the method to converge.
A value of \omega = 0.5 is a reasonable starting point.
With some trial and error you should be able to get the calculation to converge about twice as fast as
the simple relaxation method, i.e., in about half as many iterations.

d) Are there any circumstances under which using a value \omega < 0 would help us find a solution faster than
we can with the ordinary relaxation method? (Hint: The answer is yes, but why?)
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <string.h>     //string routines
#include <math.h>       //math routines

#include "relaxation.h" //relaxation method header file

void main()
{
    float f(float x);                //right hand side of the equation x = f(x)
    float relaxation(float x_0, float (*func_ptr)(float x), float accuracy);  //standard relaxation method

    float x_0;                       //starting point
    float x;
    float (*func_ptr)(float x);
    float accuracy;                  //target accuracy

    func_ptr = f;
    strcpy(equation_name, "x=1-exp(-2x)");
    printf("Solution to %s using relaxation method\n\n", equation_name);
    //printf("Enter the starting point x_0 = ");              //enter manually the starting point and the target accuracy
    //scanf("%f", &x_0);
    //printf("\nEnter the target solution accuracy \\epsilon = ");
    //scanf("%f", &accuracy);
    //printf("\n");
    x_0 = 1;                                                  //predefined starting point
    accuracy = 1e-6;                                          //predefined target accuracy
    printf("x_0 = %7.6e, accuracy = %7.6e\n", x_0, accuracy);

    x = relaxation(x_0, func_ptr, accuracy);
    printf("\nx = %7.6e\n\n", x);

    printf("Solution to %s using overrelaxation method\n\n", equation_name);
    x = overrelaxation(x_0, func_ptr, accuracy);
}

float f(float x)
{
    return (1 - exp(- 2 * x));
}

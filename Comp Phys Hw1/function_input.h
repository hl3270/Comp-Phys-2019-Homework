/*
function_input.h

To choose or input one dimensional or multi dimensional functions to be used in calculation

Required header files:
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/

/*
A set of predefined iD functions FUNCTION_1_i, with a related index i
*/
#define FUNCTION_1_1(x) cos(x)
#define FUNCTION_1_2(x) exp(x)
#define FUNCTION_1_3(x) exp(-x)

/*
Expressions of FUNCTION1_i above, used in filename and inside the file
*/
#define FUNCTION_EXPRESSION_1_1 "cos(x)"
#define FUNCTION_EXPRESSION_1_2 "exp(x)"
#define FUNCTION_EXPRESSION_1_3 "exp(-x)"

int function_choice;      //index for choice within a specific set of functions predefined
char *function_1D;        //function name string to put in the filename

/*
char *function_input()
{
    char *input_f;

    printf("Input a one dimensional function of x you want to differentiate : cos(x) or exp(x)\n");
    scanf("%s", &input_f);
    printf("\n");

    return *input_f;
}
*/

/*
Choose within a set of functions by inputing a related index
*/
int choose_function()
{
    int i;

    do
    {
        printf("Choose one function to calculate. \nFor cos(x), type 1; for exp(x), type 2; for exp(-x), type 3. \n");
        scanf("%d", &i);
        printf("\n");

        switch(i)
        {
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        default:
            {
                printf("Wrong choice! Please type 1 or 2. \n");
                i = 0;
            }
        }
    } while (i == 0);

    return i;
}

/*
Predefined one dimensional functions. for index i chosen above, here return FUNCTION1_i defined at the very above
*/
float f_1D(float x)
{
    switch(function_choice)
    {
    case 1:
        return FUNCTION_1_1(x);
        break;
    case 2:
        return FUNCTION_1_2(x);
        break;
    case 3:
        return FUNCTION_1_3(x);
    default:
        break;
    }
}

/*
To put the expression of functions into the filename and the file    ////try to input the functions manually
*/
char *f_1D_expression(int function_choice)
{
    switch (function_choice)
    {
    case 1:
        return FUNCTION_EXPRESSION_1_1;
        break;
    case 2:
        return FUNCTION_EXPRESSION_1_2;
        break;
    case 3:
        return FUNCTION_EXPRESSION_1_3;
        break;
    default:
        break;
    }
}

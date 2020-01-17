/*Hw_6-2.c*/

/*
Exercise 10.11 (Computational Physics, Mark Newman)
The dimer covering problem

A well studied problem in condensed matter physics is the dimer covering problem in which dimers,
meaning polymers with only two atoms, land on the surface of a solid, falling in the spaces
between the atoms on the surfaces and forming a grid like this.

No two dimers are allowed to overlap. The question is how many dimers we can fit in the entire L * L square.
The answer, in this simple case, is clearly L * L / 2, but suppose we did not know about this.
(There are more complicated versions of the problem on different lattices, or with differently shaped elements,
for which the best solution is far from obvious, or in some cases not known at all.)

a) Write a program to solve the problem using simulated annealing on a 50 * 50 lattice.
The "energy" function for this system is minus the number of dimers, so that is minimized
when the dimers are a maximum. The moves for the Markov chain are as follows:
  i) Choose two adjacent sites on the lattice at random
  ii) If those two sites are currently both empty, place a dimer on them.
  iii) If the two sites are currently occupied by a single dimer, remove the dimer from the lattice
       with the appropriate probability (which you will have to work out).
  iv) Otherwise, do nothing.
Create an animation of the state of the system over time as the simulator runs.

b) Try exponential cooling schedules with different time constants.
A reasonable first value to try is \tau = 10,000 steps. For faster cooling schedules you should see that
the solutions found are poorer - a smaller fraction of the lattice is filled with dimers and
there are larger holes in between them - but for slower schedules the calculation can find quite good,
but usually not perfect, coverings of the lattice.

One thing that is not clear from the instructions for the Markov chain is that,
when the pair of cells are empty you always add a dimer to them because this lowers the total energy of the system.
If the pair of cells has one dimer connecting them, you remove that dimer with probability equal to
the Boltzmann probability. Since the total energy is minus the sum of all dimers,
the probability of removing one dimer is P_acc = exp(- 1 / T),
where T is the current temperature of the system at that time in the cooling schedule.
NB, once again, you do not need to make a movie, but you should include:
a) a plot comparing an early state in the system with the final state in the system.
b) a table (or equivalent) showing how the solution depends on the cooling schedule
you have chosen.
*/

#include <stdlib.h>     //standard library
#include <stdio.h>      //standard input/output library
#include <stdbool.h>    //standard boolean library
#include <string.h>     //string routines
#include <math.h>       //math routines
#include <malloc.h>     //memory allocation routines
#include <time.h>       //time routines

#define SEED time(NULL)

double kb = 1;
double T_0 = 10;
double tau = 1000;
int N = 100000;

void main()
{
    int di[50][50];       //50 * 50 lattice showing if a site is occupied by a dimer. di = 0 if unoccupied, and 1 if occupied
    int ne[50][50];       //50 * 50 lattice showing which neighbor of this site is the other atom if it is occupied by a dimer.
                          //- 1 if unoccupied, 0 for x + 1, 1 for y + 1, 2 for x - 1, 3 for y - 1
    int rx;               //choose a random x
    int ry;               //choose a random y
    int rn;               //choose a random neighbor
    int nx;               //x of the random neighbor
    int ny;               //y of the random neighbor
    double r;             //random number
    double *E;            //total energy
    //double T_0;           //initial temperature
    //double tau;           //time constant
    double T;             //temperature
    //double kb;            //Boltzmann constant
    //int N;                //target Monte Carlo steps
    int n;                //current step
    int i;                //x index
    int j;                //y index
    int d;                //number of sites occupied
    int exit;
    FILE *fp;
    char filename[64];
    FILE *fp1;
    char filename1[64];

    //kb = 1;
    //T_0 = 10;
    //tau = 5000;
    //N = 60000;

    srand(SEED);

    /*Initialize dimer lattice*/
    for(i = 0; i < 50; i++)
    {
        for(j = 0; j < 50; j++)
            {
                di[i][j] = 0;
                ne[i][j] = - 1;
            }
    }
    E = (double *) malloc((N + 1) * sizeof(double));
    E[0] = 0;

    /*Add or remove a dimer until target step*/
    for(n = 1; n <= N; n++)
    {
        /*Get a random adjacent pair of sites*/
        rx = (int) (50 * rand() / (RAND_MAX + 1));  //get a random number from 0 to 49
        ry = (int) (50 * rand() / (RAND_MAX + 1));
        nx = rx;
        ny = ry;
        do
        {
            exit = 0;
            rn = (int) (4 * rand() / (RAND_MAX + 1));
            //printf("%d\n", rn);
            if(rn == 0)                            //if the neighbor chosen is outside the lattice, choose neighbor again
            {
                nx += 1;
                if(nx > 49)
                    exit = 1;
            }
            else if(rn == 1)
            {
                ny += 1;
                if(ny > 49)
                    exit = 1;
            }
            else if(rn == 2)
            {
                nx -= 1;
                if(nx < 0)
                    exit = 1;
            }
            else if(rn == 3)
            {
                ny -= 1;
                if(ny < 0)
                    exit = 1;
            }
        } while(exit == 1);
        //printf("%d\n", rx);

        if((di[rx][ry] == 0) && (di[nx][ny] == 0)) //if the pair of sites chosen are both unoccupied, add a dimer here
        {
            di[rx][ry] = 1;                        //add dimer
            di[nx][ny] = 1;
            ne[rx][ry] = rn;                       //label the neighbor
            if(rn == 0)
            {
                ne[nx][ny] = 2;
            }
            else if(rn == 1)
            {
                ne[nx][ny] = 3;
            }
            else if(rn == 2)
            {
                ne[nx][ny] = 0;
            }
            else if(rn == 3)
            {
                ne[nx][ny] = 1;
            }
        }
        else if(ne[rx][ry] == rn)                  //if the pair of cites are occupied by a single dimer, remove it with Boltzmann probability
        {
            r = (double) (rand() / (RAND_MAX + 1));
            T = T_0 * exp(- n / tau);
            if(r < exp(- 1 / (kb * T)))            //move accepted with Boltzmann probability
            {
                di[rx][ry] = 0;
                di[nx][ny] = 0;
                ne[rx][ry] = - 1;
                ne[nx][ny] = - 1;
            }
        }

        d = 0;
        for(i = 0; i < 50; i++)
        {
            for(j = 0; j < 50; j++)
            {
                d += di[i][j];
            }
        }
        E[n] = (double) (- d / 2);

        printf("%f\n", E[n]);

        if((n == 100) || (n == N))
        {
            sprintf(filename1, "dimer_lat_T0_%2.1e_tau_%3.1e_n_%d.txt", T_0, tau, n);
            fp1 = fopen(filename1, "w");
            fprintf(fp1, "Dimer covering lattice\n\n", T, n);
            fprintf(fp1, "initial temperature T_0 = %2.1e, time constant tau = %3.1e, n = %d\n\n", T_0, tau, n);
            for(i = 0; i < 50; i++)
            {
                for(j = 0; j < 50; j++)
                {
                    fprintf(fp1, "%d ", di[i][j]);
                }
                fprintf(fp1, "\n");
            }
            fclose(fp1);
        }
    }
    printf("E  %f %d\n", E[N], N);

    sprintf(filename, "dimer_T0_%2.1e_tau_%3.1e.txt", T_0, tau);
    fp = fopen(filename, "w");
    fprintf(fp, "Dimer covering\n\n", T, n);
    fprintf(fp, "initial temperature T_0 = %2.1e, time constant tau = %3.1e\n\n", T_0, tau);
    fprintf(fp, "E\n");
    for(int k = 0; k <= N; k++)
        fprintf(fp, "%5.4e\n", E[k]);
    fclose(fp);

    free(E);
}

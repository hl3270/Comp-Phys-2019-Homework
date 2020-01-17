/*
SBHB.h

Header file for equation of motion of a supermassive black hole binary (SBHB)
*/

double A;
double B;
int dim = 4;

double r_s = 1e-7;

double *f(int d, double *x, double t)
{
    double *y;
    double r;
    double v;

    y = (double *) malloc(d * sizeof(double));

    r = sqrt(x[0] * x[0] + x[1] * x[1]);
    v = sqrt(x[2] * x[2] + x[3] * x[3]);

    y[0] = x[2];
    y[1] = x[3];
    y[2] = - x[0] / (4 * r * r * r) - A * x[2] / (B + v * v * v);
    y[3] = - x[1] / (4 * r * r * r) - A * x[3] / (B + v * v * v);

    return y;
}

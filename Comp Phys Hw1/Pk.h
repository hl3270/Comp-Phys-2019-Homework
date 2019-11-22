/*
Pk.h

Functions and variables specially for power spectrum P(k) data and data files
*/

#include "data_file.h"         //data and file related header file
#include "cubic_spline.h"      //cubic spline header file
#include "integral_1D.h"       //iD integral header file

#define DATA_FILE_SUFFIX ".matter_pk"                   //suffix of filename of the data file type we are dealing with in this problem

#define P_1K(K) (Pk[0].P * pow((K / Pk[0].k), dfl))     //P(k) by estimation if cubic spline, estimate log(P)-log(k) to be linear out of boundary
#define P_2K(X) (Pk[J].P * pow(10, ((cs_logPk.a[J] * pow((X - logPk.x[J]), 3)) + (cs_logPk.b[J] * pow((X - logPk.x[J]), 2)) + cs_logPk.c[J] * (X - logPk.x[J]))))  //X = log10(k)
//#define P_2K(X) (Pk[J].P)
//#define P_3K(Y) (Pk[N - 1].P * pow((1 / (Y * Pk[N - 1].k)), dfr))  //Y = 1/k, P(k) by estimation if cubic spline, estimate log(P)-log(k) to be linear out of boundary
#define F_PK_1(K) ((K * P_1K(K) * sin(K * r)) / (2 * pow(PI, 2) * r))
#define F_PK_2(X) (log(10) * pow(100, (X)) * P_2K(X) * sin(r * pow(10, X)) / (2 * pow(PI, 2) * r))
//#define F_PK_3(Y) ((- 1 / pow(Y, 3)) * P_3K(Y) * sin(r / Y) / (2 * pow(PI, 2) * r))
#define F_PK_3(K) (K * Pk[N - 1].P * pow((K / Pk[N - 1].k), dfr) * sin(r * K) / (2 * pow(PI, 2) * r))

#define PI (3.14159265)        //value of \pi

int N;                         //the number of rows of data in the file, i.e. the number of data points

struct two_dim_data logPk;     //log(k)-log(P) data for cubic spline to estimate P(k) for all k
struct cspline cs_logPk;       //natural cubic spline result of log(k)-log(P)
float dfl;                     //derivative of log(P)-log(k) at k_0
float dfr;                     //derivative of log(P)-log(k) at k_n, n = N - 1
float r;                       //scale r
int J;                         //index of bins of log(P)-log(k)
int fPk;                       //for a specific r, return a related function depend on range of k, for i, return F_PK_i

/*
structure of P(k) data stored in ".matter_pk" files and to be calculated
*/
struct matter_pk
{
    float k;                   //wavenumber
    float P;                   //power spectral density
    float third_quantity;      //unknown quantity on the third column
    float fourth_quantity;     //unknown quantity on the fourth column
};

struct matter_pk *Pk;          //data to deal with

/*
structure of \xi(r) calculated from P(k)
*/
struct xi_r_pk
{
    int n;                     //number of r's
    float *r;                  //sequence of r
    float *xi;                 //sequence of \xi of related r
};

struct xi_r_pk xi_r_Pk;        //result of correlation function

/*
Input power spectrum P(k) manually

Structure of data:
wavenumber k, power spectral density P, unknown quantity on the third column, and unknown quantity on the fourth column
*/
struct matter_pk *input_manual_Pk() ///////////////////unfinished
{
    struct matter_pk *pk;

    return (pk);
}

/*
Input power spectrum data P(k) by file

Structure of data:
wavenumber k, power spectral density P, unknown quantity on the third column, and unknown quantity on the fourth column
*/
struct matter_pk *input_file_Pk(char *filename, int N)
{
    struct matter_pk *pk;
    FILE *fp;

    /*Read the data in the file and transfer as input*/
    pk = (struct matter_pk *)malloc(N * sizeof(struct matter_pk));
    fp = fopen(filename, "r");
    for (int i = 0; i < N; i++)      //read P(k) data in structure row by row from data file.
    {
        fscanf(fp, "%e %e %e %e\n", &pk[i].k, &pk[i].P, &pk[i].third_quantity, &pk[i].fourth_quantity);
    }
    fclose(fp);

    return pk;
}

/*
Input power spectrum data P(k), either by file or manually
*/
void input_Pk()
{
    int Data_Input_Method();         //choose a method to input data. 1: by file; 2: by hand
    struct FileName open_data_file(char *data_file_suffix);   //open a file of the specific type with the filename entered and save the filename
    int row_file(char *filename);    //count the number of rows of data in the file opened
    struct matter_pk *input_file_Pk(char *filename, int N);   //input data from the ".matter_pk" file opened
    struct matter_pk *input_manual_Pk();                      //input Pk data manually

    int data_input_method;           //method to input data: 1: by file; 2: by hand

    data_input_method = Data_Input_Method();                  //choose method to input data. 1 : input by file; 2 : input manually
    if (data_input_method == 1)
    {
        fn = open_data_file(DATA_FILE_SUFFIX);
        //printf("%s", fn.filename);
        N = row_file(fn.filename);
        //printf("%d\n",N);
        Pk = input_file_Pk(fn.filename, N);
        //printf("%e", Pk[0].P);
    }
    else
        Pk = input_manual_Pk();

    //for(i = 0; i <= N-1; i++) printf("%7.6e\n", Pk[i].P);
}

/*
Calculate correlation function \xi(r) with P(k) data
*/
void corr_func(struct matter_pk *Pk, struct cspline cs_logPk)
{
    float Simpson(float a, float b, float (*func_ptr)(float), int N_0, float esp_tol, int max_loop); //Simpson's method integration
    float Romberg(float a, float b, float (*func_ptr)(float), int N_0, float esp_tol, int max_loop); //Romberg integration
    float f(float x);              //function to integrate when calculating correlation function

    float r_lower;             //lower bound of the scale to calculate
    float r_upper;             //upper bound of the scale to calculate
    int n_r;                   //bins of r, i.e. divide the range of r to n_r identical parts
    int k_upper;               //user defined upper bound of wavenumber k. For infinite integral of highly oscillatory functions, it is not convenient even if we change the variable to, e.g. 1 / k
    float xi_1;                //correlation integral over [0, k_0]
    float xi_2;                //correlation integral over [k_0, k_n]
    float xi_3;                //correlation integral over [k_n, +INF], since the function to calculate is highly oscillatory no matter change the variable or not, we simply integrate to a user defined upper bound of k and ignore the tail
    float xi;                  //correlation function value
    //struct xi_r_pk xi_r;       //correlation function result

    char filename[64];         //filename string
    FILE *fp;                  //pointer to file to store

    dfl = cs_logPk.c[0];
    dfr = 3 * cs_logPk.a[cs_logPk.n - 1] * pow((cs_logPk.x[cs_logPk.n] - cs_logPk.x[cs_logPk.n - 1]), 2) + 2 * cs_logPk.b[cs_logPk.n - 1] * (cs_logPk.x[cs_logPk.n] - cs_logPk.x[cs_logPk.n - 1]) + cs_logPk.c[cs_logPk.n - 1];
    //printf("%7.6e %7.6e\n", dfl, dfr);

    /*Enter parameters of the calculation of correlation function*/
    printf("Enter the range of scale r you want to calculate : [r_lower, r_upper] \nr_lower = ");
    scanf("%f", &r_lower);
    printf("r_upper = ");
    scanf("%f", &r_upper);
    printf("How many parts would you like to separate the range of r in calculation : n_r = ");
    scanf("%d", &n_r);
    printf("Enter for ease of calculation an upper bound of wavenumber k : ");
    scanf("%d", &k_upper);

    sprintf(filename, "corr_func_from_%s_within_[%f,%f].txt", fn.filename, r_lower, r_upper);
    fp = fopen(filename, "w");
    fprintf(fp, "Correlation function xi(r) calculated from %s within r = [%f, %f]\n\n", fn.filename, r_lower, r_upper);
    fprintf(fp, "r : scale; xi_1 : integral between 0 and k_min; xi_2 : integral within data range; xi_3 : integral from k_max to +INF\n\n");
    fprintf(fp, "    r        xi_1         xi_2         xi_3          xi       r^2*xi_2      r^2*xi\n");

    //xi_r.n - n_r + 1;
    for (int i = 0; i <= n_r; i++)
    {
        r = r_lower + i * (r_upper - r_lower) / n_r;

        fPk = 1;               //integrate over k \in [0, k_min], with the assumed P(k) from the lower boundary of cubic spline result
        xi_1 = Romberg(0, Pk[0].k, f, 2, 1e-7, 5);
        //printf("%7.6e\n", xi_1);

        fPk = 2;               //integrate over the data range with the cubic spline result
        xi_2 = 0;
        for(J = 0; J <= logPk.n - 1; J++)
        {
            xi_2 += Simpson(logPk.x[J], logPk.x[J + 1], f, 2, 1e-7, 12);
        }
        //printf("%7.6e\n", xi_2*r*r);

        fPk = 3;               //integrate over k \in [k_max, +INF), with the assumed P(k) from the upper boundary of the cubic spline result
        xi_3 = Simpson(Pk[N - 1].k, k_upper, f, 500, 1e-7, 12);
        //printf("%7.6e\n", xi_3);

        //xi_r.r[i] = r;
        //xi_r.xi[i] = xi_1 + xi_2 + xi_3;

        xi = xi_1 + xi_2 + xi_3;

        fprintf(fp, "%6.2f %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n", r, xi_1, xi_2, xi_3, xi, r * r * xi_2, r * r * xi);
        printf("%7.6e\n", r * r * xi);

        r += (r_upper - r_lower) / n_r;
    }

    fclose(fp);
};

/*
Functions to integrate for different range of wavenumber k
*/
float f(float x)
{
    switch (fPk)
    {
    case 1:
        return F_PK_1(x);      //k \in [0, k_min], P(k) \propto k ^ a, a > 0
    case 2:
        return F_PK_2(x);      //k \in [k_min, k_max], log(P)-log(k) follows the cubic spline result in between k points
    case 3:
        return F_PK_3(x);      //k \in [k_max, +INF], P(k) \propto k ^ a, a < 0
    default:
        break;
    }
}

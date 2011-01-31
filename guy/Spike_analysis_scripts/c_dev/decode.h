/* Header file for decode function
 * Guy Billings UCL 2010
 */
 
#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define __GSL_RANGE_CHECK_OFF__
extern gsl_rng *rgen;

void import_gsl_matrix(gsl_matrix *,double *,int,int);
double custom_rand(double,gsl_rng *);
double * build_double_array(int);
void import_gsl_matrix(gsl_matrix *,double *,int,int);
double custom_rand(double,gsl_rng *);
double db_euclid(int,int,int,gsl_matrix *,gsl_matrix *);
double find_closest(int, double *, gsl_rng *);
double retrieve_value(int, int, int, double *);
gsl_matrix* build_gsl_matrix(int, int, double, int);
void copy_double_arr(double *, double *, int);
 
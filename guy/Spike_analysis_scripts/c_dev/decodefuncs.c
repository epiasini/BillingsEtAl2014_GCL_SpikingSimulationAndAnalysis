/* Functions used in the decode function
 * Guy Billings UCL, 2010
 */
 
#include "decode.h"

/* Function to allow the creation of dynamic 1d array (double) */
double * build_double_array(int done)
{
    int i;
    double *arr;
    arr=malloc(done * sizeof(double));
    if(arr == NULL)
    {
     /* Memory full */
     mexPrintf("Memory error!");
    }
                 
    for(i = 0; i < done ; i++) {
    arr[i] = 0.0;
    }
                
 return arr;

}

/* Matrix creation function using GSL */
gsl_matrix* build_gsl_matrix(int done, int dtwo, double init, int mode)
{
 int i,j;
 extern gsl_rng *rgen;
 gsl_matrix* m;
 m = gsl_matrix_alloc (done, dtwo);
 if(m==NULL){
   error_out(2);
 }
 else{  
 if(mode==0){ /*For mode = 0 fill with init*/
  for (i = 0; i < done; i++)
    for (j = 0; j < dtwo; j++)
      gsl_matrix_set (m, i, j,init);
  }
 else {     
 for (i = 0; i < done; i++)
    for (j = 0; j < dtwo; j++)
      gsl_matrix_set (m, i, j,custom_rand(init,rgen));
  }
 }
 return m;     
}                 
 
/* Function resolve data pointers to the matlab workspace in an 
 *easier matrix style form */
double retrieve_value(int i, int j, int done, double *ptr)
{
 double value;
 int addr=(j-1)*done+(i-1);
 value=*(ptr+addr);
 return value;
}

 /* Function to copy data from MATLAB matricies */
 void import_gsl_matrix(gsl_matrix *out, double *in, int k, int l)
 {
  unsigned register i, j;
  
  for(i=0 ; i<k ; i++)
  {
   for(j=0 ; j<l ; j++)
   {
    gsl_matrix_set(out,i,j,retrieve_value(i+1,j+1,k,in));
   }
  }
 }   
  
/* Function that returns a random value between two bounds drawn from a      */
/* uniform distribution using the currently selected random number generator */
double custom_rand(double maxi, gsl_rng *generator)
{
 double varate;
 /* use high order bits */
 varate=maxi*gsl_rng_get(generator)/gsl_rng_max(generator);
 /*varate=maxi*(rand()/(RAND_MAX+1.0));*/
 return varate;
}

/* Function to calculate the Euclidean distance an observation and a cluster */
double db_euclid(int c,int o,int vector_len,gsl_matrix *book_vector,gsl_matrix *data_vector)
{
 unsigned register i;
 double sumsquare=0;
 
 for(i=0 ; i<vector_len ; i++)
 {  
    sumsquare=sumsquare+pow(gsl_matrix_get(book_vector,c,i)-gsl_matrix_get(data_vector,o,i),2);
 }
  
  return sqrt(sumsquare);
}

/* Function to find the closest cluster. If more than one is at */
/* equal distance, then a random one is chosen.					*/
double find_closest(int clusters, double *dist, gsl_rng *generator)
{
  unsigned register i;

  double min=*dist;
  double minc=0;

  for (i=1; i<clusters; i++)
  {
   if (*(dist+i)<min)
     min=(double)*(dist+i);
     minc=(double)i;
  }
  return minc;
}  
      
/* Function to copy arrays in to matlab ones */
 void copy_double_arr(double *in, double *out, int size)
 {
  int i;
  
  for(i=0 ; i<size ; i++)
  {
    *(out+i)=in[i];
  }
 }      
      
   
   
      



 
     
 


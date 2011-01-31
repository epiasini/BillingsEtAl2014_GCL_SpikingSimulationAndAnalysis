/* MEX decoder function
 * Guy Billings UCL 2010
 * Usage: [alphabet]=decode(observations,clusters,book_vector,data_vector,vector_len,debug)
 * alphabet : 		Array of values with index:observation and value:cluster 
 * observations:	Integer number of data points
 * clusters:		Integer number of clusters
 * book_vector:		Code book vector for each cluster, rows:clusters columns:data
 * data_vector:		Vector for each observation, rows:observations columns:data
 * vector_len:		Length of each data/book vector
 */
#include "decode.h"

gsl_rng *rgen; /* random number generator is global */

void mexFunction(int nlhs,      mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    if (nrhs!=6)
    {
     mexPrintf("Wrong number of arguments to decodeFST\n");
     return;
    } 
    
    double debug=mxGetScalar(prhs[5]);

    if(debug>=1)
      mexPrintf("Intialising: Random number generator\n");
      
	/*initialise random number generator*/
	rgen=gsl_rng_alloc(gsl_rng_taus);
	
	if(debug>=1)
      mexPrintf("Intialising: Interface variables\n");
	
    /*interface variables */
    const int observations=(int)mxGetScalar(prhs[0]);
    const int clusters=(int)mxGetScalar(prhs[1]);
    const int vector_len=(int)mxGetScalar(prhs[4]);
    
    if(debug>=1)
      mexPrintf("Intialising: GSL matricies\n");
    
    /*use gsl matricies for input matrix data */
    gsl_matrix *book_vector;
    book_vector=build_gsl_matrix(clusters,vector_len,0,0);
    gsl_matrix *data_vector;
    data_vector=build_gsl_matrix(observations,vector_len,0,0);
    
    if(debug>=1)
      mexPrintf("Intialising: Output variable\n");
    
    /*initialise output variable */
    double *alphabet;
    alphabet=build_double_array(observations);
    double *alph;
    
    if(debug>=1)
      mexPrintf("Intialising: Internal data\n");
    
    /*initialise internal variables */
    double *dist;
    dist=build_double_array(clusters);
    unsigned register o,c;
    
    if(debug>=1)
      mexPrintf("Loading matrix data\n");
    
    /*load matrix data */
    import_gsl_matrix(book_vector,mxGetPr(prhs[2]),clusters,vector_len);
    import_gsl_matrix(data_vector,mxGetPr(prhs[3]),observations,vector_len);
    

    for(o=0 ; o<=observations-1 ; o++)
    {
      /*for each observation determine distance to each cluster center */
      for(c=0 ; c<=clusters-1 ; c++)
      { 
        *(dist+c)=db_euclid(c,o,vector_len,book_vector,data_vector);
        if(debug>=1)
          mexPrintf("Distance %5.4e from cluster %d to observation %d\n",*(dist+c),c+1,o+1);
      }
      
      /*find the nearest cluster and assign that cluster to the observation*/
      *(alphabet+o)=find_closest(clusters,dist,rgen)+1;
      if(debug>=1)
         mexPrintf("Nearest cluster to observation %d is %d\n",o+1,*(alphabet+o));
    }  
    
    if(debug>=1)
      mexPrintf("Copying back to workspace\n");
    /*copy results back to workspace */
    plhs[0]=mxCreateDoubleMatrix(1,observations,mxREAL);
    alph=mxGetPr(plhs[0]);
    copy_double_arr(alphabet,alph,observations);
    if(debug>=1)
      mexPrintf("Clearing up memory\n");
    /*deallocate memory */
    free(dist);
    free(alphabet);
    gsl_matrix_free(book_vector);
    gsl_matrix_free(data_vector);
  
}    
    
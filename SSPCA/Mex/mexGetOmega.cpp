#include "mex.h"
#include <cmath>

  // mex mexGetOmega.cpp -largeArrayDims

typedef ptrdiff_t intt;

  // Argument 0: V,    whose dim is p x r
  //          1: Gsq,  whose dim is p x ng (SPARSE)
  //          2: params

//-------------------------------------------------------------------------

void parse_params(const mxArray *params, double *normparam);

//-------------------------------------------------------------------------

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
  
  double *V, *Gsq, *Omega;
  intt   *G_ir, *G_jc;;

  int p = mxGetM(prhs[0]), r = mxGetN(prhs[0]), ng = mxGetN(prhs[1]);

  V     = (double*)mxGetPr(prhs[0]);
  Gsq   = (double*)mxGetPr(prhs[1]);
 
  G_ir  = (intt*)mxGetIr(prhs[1]);
  G_jc  = (intt*)mxGetJc(prhs[1]);
  
  
  plhs[0]  = mxCreateDoubleScalar( 0.0 ); //mxCreateDoubleMatrix(1, 1, mxREAL);//
  Omega    = mxGetPr( plhs[0] );
  //-----------------------------------------------------------------------
  double normparam;
  
  parse_params(prhs[2], &normparam);
  
  const double power1 = normparam / 2.0, power2 = 1 / normparam;
  //-----------------------------------------------------------------------
  double *V_sq          = new double[p*r];

  for (intt i = 0; i<p*r; i++)
    V_sq[i] = V[i] * V[i];  
  //-----------------------------------------------------------------------
  double block, norm_for_one_col_of_V;
  
  for (intt col = 0; col < r; col++)
    {
       
    norm_for_one_col_of_V = 0.0;
    
    for (intt g = 0; g < ng; g++)
      {
          
       block = 0.0; 
        
       for (intt idx = G_jc[g]; idx < G_jc[g+1]; idx++)
            block += V_sq[col*p+G_ir[idx]] * Gsq[idx];
      

       norm_for_one_col_of_V += pow( block, power1 );

      }

     *Omega += pow( norm_for_one_col_of_V, power2 );;

    }
  //-----------------------------------------------------------------------
  delete [] V_sq;
}

//=========================================================================

void parse_params(const mxArray *params, double *normparam)
{
   //-----------------------------------------------------------------------
   mxArray *tmp_field; 
   mwIndex index_params = 0;

   //-----------------------------------------------------------------------
   tmp_field = mxGetField(params, index_params, "normparam");

   if (tmp_field == NULL)
      *normparam = 0.5;
   else
      *normparam = (double)*mxGetPr(tmp_field); 
   //-----------------------------------------------------------------------  
}

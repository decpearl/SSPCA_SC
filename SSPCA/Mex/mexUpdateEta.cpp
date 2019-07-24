#include "mex.h"
#include <cmath>

  // mex mexUpdateEta.cpp -largeArrayDims

typedef ptrdiff_t intt;

  // Argument 0: V,    whose dim is p x r
  //          1: Gsq,  whose dim is p x ng (SPARSE)
  //          2: params

//-------------------------------------------------------------------------

void parse_params(const mxArray *params, double *epsilon, double *normparam);

//-------------------------------------------------------------------------

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
  
  double *V, *Eta, *inv_Zeta, *Gsq;
  intt   *G_ir, *G_jc;;

  int p = mxGetM(prhs[0]), r = mxGetN(prhs[0]), ng = mxGetN(prhs[1]);

  V     = (double*)mxGetPr(prhs[0]);
  Gsq   = (double*)mxGetPr(prhs[1]);
 
  G_ir  = (intt*)mxGetIr(prhs[1]);
  G_jc  = (intt*)mxGetJc(prhs[1]);
  
  //-----------------------------------------------------------------------
  double epsilon, normparam;
  
  parse_params(prhs[2], &epsilon, &normparam);
  
  const double power1 = (2.0-normparam)/2.0, power2 = normparam/2.0, power3 = (normparam-1.0)/normparam;
  //-----------------------------------------------------------------------
  plhs[0]  = mxCreateDoubleMatrix(p, r, mxREAL);
  inv_Zeta = (double*)mxGetPr( plhs[0] );

  plhs[1]  = mxCreateDoubleMatrix(ng, r, mxREAL);
  Eta      = (double*)mxGetPr( plhs[1] );
  
  for (intt i = 0; i < p*r; i++)
    inv_Zeta[i] = 0.0;
  
  for (intt i = 0; i < ng*r; i++)
    Eta[i] = 0.0;
  //-----------------------------------------------------------------------
  double *V_sq          = new double[p*r];
  double *Normalization = new double[r];

  for (intt i = 0; i<p*r; i++)
    V_sq[i] = V[i] * V[i];  

  for (intt i = 0; i < r; i++)
    Normalization[i] = 0.0;
  //-----------------------------------------------------------------------
  for (intt col = 0; col < r; col++)
    {
    for (intt g = 0; g < ng; g++)
      {
          
       for (intt idx = G_jc[g]; idx < G_jc[g+1]; idx++)
            Eta[col*ng+g] += V_sq[col*p+G_ir[idx]] * Gsq[idx];
      

       Normalization[col] += pow( Eta[col*ng+g], power2 );

       Eta[col*ng+g]       = pow( Eta[col*ng+g], power1 );

      }

      if ( Normalization[col] != 0 )
        Normalization[col] = pow( Normalization[col], power3 );

    }
  //-----------------------------------------------------------------------

  double inv_eta;

  for (intt col = 0; col < r; col++)
    for (intt g = 0; g < ng; g++)
      {

        Eta[col*ng+g] *= Normalization[col];//Eta[col*ng+g] += epsilon;
        
        if ( Eta[col*ng+g] < epsilon )
            Eta[col*ng+g] = epsilon;
        
        inv_eta = 1.0 / Eta[col*ng+g];
	
        
        for (intt idx = G_jc[g]; idx < G_jc[g+1]; idx++)
            inv_Zeta[col*p+G_ir[idx]] += Gsq[idx] * inv_eta;        
        
      }


  delete [] Normalization;
  delete [] V_sq;
}

//=========================================================================

void parse_params(const mxArray *params, double *epsilon, double *normparam)
{
   //-----------------------------------------------------------------------
   mxArray *tmp_field; 
   mwIndex index_params = 0;

   //-----------------------------------------------------------------------
   tmp_field = mxGetField(params, index_params, "epsilon"); 

   if (tmp_field == NULL)
      *epsilon = 1e-5;
   else 
      *epsilon = (double)*mxGetPr(tmp_field);
   //-----------------------------------------------------------------------
   tmp_field = mxGetField(params, index_params, "normparam");

   if (tmp_field == NULL)
      *normparam = 0.5;
   else
      *normparam = (double)*mxGetPr(tmp_field); 
   //-----------------------------------------------------------------------  
}

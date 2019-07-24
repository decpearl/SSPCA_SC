#include "mex.h"
#include "blas.h"

// mex mexUpdateV.cpp -l blas -largeArrayDims

typedef ptrdiff_t intt;

// Argument 0: V, whose dim is p x r
//          1: X, whose dim is n x p
//          2: U, whose dim is n x r
//          3: inv_Zeta, whose dim is p x r
//          4: lambda
//          5: params

//-------------------------------------------------------------------------

void parse_params(const mxArray *params, intt *max_it, bool *pos, double *normalization);

//=========================================================================

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

   double *U, *V, *X, *inv_Zeta, *V_output, lambda;

   double one_db = 1.0, zero_db = 0.0, minus_one_db = -1.0;

   intt one_intt = 1, zero_intt = 0;
   intt p  = mxGetM(prhs[0]), r =  mxGetN(prhs[0]), n = mxGetM(prhs[1]);
   intt pr = p*r;

   V        = (double*)mxGetPr(prhs[0]);
   X        = (double*)mxGetPr(prhs[1]);
   U        = (double*)mxGetPr(prhs[2]);
   inv_Zeta = (double*)mxGetPr(prhs[3]);
   lambda   = (double)*mxGetPr(prhs[4]);
   //-----------------------------------------------------------------------
   intt   max_it;
   bool   pos;
   double normalization;

   parse_params( prhs[5], &max_it, &pos, &normalization);
   //-----------------------------------------------------------------------
   //V_output that is initialized as a copy of V
   plhs[0]  = mxCreateDoubleMatrix(p, r, mxREAL);
   V_output = (double*)mxGetPr( plhs[0] );

   dcopy_ ( &pr, V, &one_intt, V_output, &one_intt );
   //-----------------------------------------------------------------------
   double *XtU = new double[p*r];
   double *UtU = new double[r*r];

   for (intt i; i < p*r; i++ )
      XtU[i] = 0.0;
   for (intt i; i < r*r; i++ )
      UtU[i] = 0.0;

   // XtU <- X^T*U
   dgemm_ ("T","N",&p,&r,&n, &one_db,X,&n,U,&n,&zero_db, XtU,&p);
 
   //UtU <- U^T*U
   dsyrk_ ( "L", "T", &r, &n, &one_db, U, &n, &zero_db, UtU, &r );
   
   for (intt i = 0; i < r; i++)
       for (intt j = i+1; j < r; j++)
            UtU[j*r+i] = UtU[i*r+j];
   //-----------------------------------------------------------------------
   double *tmp = new double[p];
   double UtU_kk;
   
   
   for (intt t = 1; t <= max_it; t ++ )
      for (intt col = 0; col < r; col++) 
      {
       
       //------------------------------------------------------------------
       
        UtU_kk = UtU[col + r*col];
        
        if ( UtU_kk > 1e-12 )
        {
       
        // Tmp <-   XtU(:,k)
        dcopy_ ( &p, XtU+col*p, &one_intt, tmp, &one_intt );
        
        // Tmp <-   -V * UtU(:,k) + Tmp
        dgemv_ ("N", &p, &r, &minus_one_db, V_output, &p, UtU+col*r, &one_intt, &one_db, tmp, &one_intt );
       
        // Tmp <- Tmp + UtU_kk*V_output(:,k);
        daxpy_ (&p, &UtU_kk, V_output+col*p, &one_intt, tmp, &one_intt );  
        
        //-----------------------------------------------------------------
        if ( pos ) 
        {  
            for ( intt j=0; j < p; j++ )
                if ( tmp[j] < 0 )
                    V_output[col*p+j] = 0.0;
                else
                    V_output[col*p+j] = tmp[j]/( lambda*normalization*inv_Zeta[col*p+j] + UtU_kk );

        }
        else
            for ( intt j=0; j < p; j++ )    
                V_output[col*p+j] = tmp[j]/( lambda*normalization*inv_Zeta[col*p+j] + UtU_kk);
        //-----------------------------------------------------------------
        }
      }

   delete [] tmp;
   delete [] XtU;
   delete [] UtU;

}

//=========================================================================

void parse_params(const mxArray *params, intt *max_it, bool *pos, double *normalization)
{
   //-----------------------------------------------------------------------
   mxArray *tmp_field; 
   mwIndex index_params = 0;

   //-----------------------------------------------------------------------
   tmp_field = mxGetField(params, index_params, "max_it"); 

   if (tmp_field == NULL)
      *max_it = 5;
   else 
      *max_it = (intt)*mxGetPr(tmp_field);
   //-----------------------------------------------------------------------
   tmp_field = mxGetField(params, index_params, "pos");

   if (tmp_field == NULL)
      *pos = false;
   else
      *pos = (bool)*mxGetPr(tmp_field); 
   //----------------------------------------------------------------------- 
   tmp_field = mxGetField(params, index_params, "normalization");

   if (tmp_field == NULL)
      *normalization = 1.0;
   else
      *normalization = (double)*mxGetPr(tmp_field); 
   //-----------------------------------------------------------------------
}
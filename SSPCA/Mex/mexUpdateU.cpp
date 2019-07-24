#include "mex.h"
#include <cmath>
#include "blas.h"

// mex mexUpdateU.cpp -l blas -largeArrayDims

typedef ptrdiff_t intt;

// for k=1:r, U^k <- Proj_{\ell_2}( U^k + \|V^k\|_2^2 ( XV^k - UVtV^k ) )

// Argument 0: V, whose dim is p x r
//          1: X, whose dim is n x p
//          2: U, whose dim is n x r
//          3: params

//-------------------------------------------------------------------------

void parse_params(const mxArray *params, intt *max_it, bool *pos);

//-------------------------------------------------------------------------

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

   double *U, *V, *X, *U_output;
   double inv_VtV_kk, inv_norm2; 

   double one_db = 1.0, zero_db = 0.0, minus_one_db = -1.0;

   intt one_intt = 1, zero_intt = 0;
   intt p  = mxGetM(prhs[0]), r =  mxGetN(prhs[0]), n = mxGetM(prhs[1]);
   intt nr = n*r;

   V  = (double*)mxGetPr(prhs[0]);
   X  = (double*)mxGetPr(prhs[1]);
   U  = (double*)mxGetPr(prhs[2]);

   //-----------------------------------------------------------------------
   intt  max_it;
   bool pos;

   parse_params( prhs[3], &max_it, &pos );
   //-----------------------------------------------------------------------
   //U_output that is initialized as a copy of U
   plhs[0]  = mxCreateDoubleMatrix(n, r, mxREAL);
   U_output = (double*)mxGetPr( plhs[0] );

   dcopy ( &nr, U, &one_intt, U_output, &one_intt );
   //-----------------------------------------------------------------------
   double *XV  = new double[n*r];
   double *VtV = new double[r*r];

   for (intt i; i < n*r; i++ )
      XV[i] = 0.0;
   for (intt i; i < r*r; i++ )
      VtV[i] = 0.0;

   // XV <- X*V
   dgemm ("N","N",&n,&r,&p, &one_db,X,&n,V,&p,&zero_db, XV,&n);

   
   //VtV <- V^T*V
   dsyrk ( "L", "T", &r, &p, &one_db, V, &p, &zero_db, VtV, &r );
   
   for (intt i = 0; i < r; i++)
       for (intt j = i+1; j < r; j++)
            VtV[j*r+i] = VtV[i*r+j];
   //-----------------------------------------------------------------------
   double *tmp = new double[n];
   
   for (intt t = 1; t <= max_it; t ++ )
      for (intt col = 0; col < r; col++) 
      {
       
         if ( VtV[col + r*col] > 1e-12 )
         {
            inv_VtV_kk  = 1.0 / VtV[col + r*col];

            // Tmp <-   X_V(:,k)
            dcopy ( &n, XV+col*n, &one_intt, tmp, &one_intt );

            // Tmp <-   -U * VtV(:,k) + Tmp
            dgemv ("N", &n, &r, &minus_one_db, U_output, &n, VtV+col*r, &one_intt, &one_db, tmp, &one_intt );
            
            // U_output(:,k) <- Tmp/VtV(k,k) + U_output(:,k);
            daxpy (&n, &inv_VtV_kk, tmp, &one_intt, U_output+col*n, &one_intt );

            //---------------------------------------------------------
            if ( pos ) 
            {  
               inv_norm2 = 0.0;
               
               // Threshold negative components and compute norm2
               for (intt i=0; i<n; i++) 
               {
                  if (U_output[col*n+i] < 0 )
                     U_output[col*n+i] = 0.0;
                  else
                     inv_norm2 += U_output[col*n+i] * U_output[col*n+i];
               }
               inv_norm2 = 1.0 / sqrt( inv_norm2 );

            }
            else
               inv_norm2 = 1.0 / dnrm2( &n, U_output+col*n, &one_intt ); 
            
            //---------------------------------------------------------

            if ( inv_norm2 < 1.0 )
               dscal( &n, &inv_norm2, U_output+col*n, &one_intt );

            //---------------------------------------------------------

         }
      }

   delete [] tmp;
   delete [] XV;
   delete [] VtV;

}

//=========================================================================

void parse_params(const mxArray *params, intt *max_it, bool *pos)
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
}

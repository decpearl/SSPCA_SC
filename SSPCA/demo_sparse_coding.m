%==========================================================================
% Example where Omega is taken to be the standard L1 norm
% with different possible variants
%==========================================================================
randn('state',0);

n = 1000;
p = 250;
r = 36;

X = randn(n,p);

echo on

clear params

params.r              = r;
params.it0            = 5;% Cost function displayed every 5 iterations
params.min_delta_cost = 0.1;% Stopping criterion: relative decrease in the cost function smaller than 0.1
params.lambda         = 2^-20;

%--------------------------------------------------------------------------
% Required for the L1 norm setting
params.normparam  = 1;
spG               = speye(p);
%--------------------------------------------------------------------------

[ U, V ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% With positivity constraints
params.posV  = 1;
params.posU  = 1;
%--------------------------------------------------------------------------

[ U_pos, V_pos ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% Sparse coding with Lq quasi-norm regularization, q in (0,1)
q = 0.5;
params.normparam  = q;
params.posV       = 0;
params.posU       = 0;
%--------------------------------------------------------------------------

[ U_q, V_q ] = sspca( X, spG, params );

echo off
%==========================================================================
% Example where Omega is taken to be the group Lasso norm (i.e., mixed L1-L2),
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
params.lambda         = 2^-12;

%--------------------------------------------------------------------------
% Required for the L1-L2 norm setting
params.normparam  = 1;
% We consider a partition in 5 groups of 50 variables
G            = zeros(p,5);
G(1:50,1)    = 1;
G(51:100,2)  = 1;
G(101:150,3) = 1;
G(151:200,4) = 1;
G(201:250,5) = 1;

spG          = sparse(G);
%--------------------------------------------------------------------------

[ U, V ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% % Group sparse coding with Lq-L2 quasi-norm regularization, q in (0,1)
% q = 0.5;
% params.normparam  = q;
% params.lambda     = 2^-15;
% %--------------------------------------------------------------------------
% 
% [ U_q, V_q ] = sspca( X, spG, params );

echo off
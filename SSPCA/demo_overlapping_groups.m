%==========================================================================
% In this demo, we consider sets of overlapping groups of variables to
% define the mixed norm Omega, as further detailed in
%
%   (2009) R. Jenatton, G. Obozinski and F. Bach. Structured sparse principal component analysis.
%   (2009) R. Jenatton, J.-Y. Audibert and F. Bach. Structured variable selection with sparsity-inducing norms. 
%==========================================================================
rand('state',0);

n = 1000;
p = 300;
r = 36;

X = rand(n,p);

echo on
%--------------------------------------------------------------------------
% We assume that the p=300 variables are organized in a 20x15 2-dimensional grid.
% We further assume that the nonzero patterns of the columns of V are ***rectangles*** on this grid:

Height = 20;
Width  = 15;

clear mode
mode.rectangle = true;

spG = get_groups( mode, [ Height, Width ] );
%--------------------------------------------------------------------------
% We consider the Lq-L2 quasi-norm regularization, q in (0,1)

clear params

q = 0.5;
params.normparam      = q;
params.r              = r;
params.it0            = 25;% Cost function displayed every 25 iterations
params.min_delta_cost = 0.1;% Stopping criterion: relative decrease in the cost function smaller than 0.1
params.lambda         = 2^-22;

[ U, V ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% For example, we display dictionary elements 1,10 and 20

imagesc( reshape( V(:,1), Height, Width )' )

pause

imagesc( reshape( V(:,10), Height, Width )' )

pause

imagesc( reshape( V(:,20), Height, Width )' )

pause

%==========================================================================
%==========================================================================
% We now assume that the nonzero patterns of V can be diamond-shaped (with
% up to 8 faces), still on the same 20x15 2-dimensional grid:

clear mode
mode.rectangle = true;
mode.pi4       = true;

spG = get_groups( mode, [ Height, Width ] );

%--------------------------------------------------------------------------
% We still consider the Lq-L2 quasi-norm regularization, q in (0,1)

clear params

q = 0.5;
params.normparam      = q;
params.r              = r;
params.it0            = 25;% Cost function displayed every 25 iterations
params.min_delta_cost = 0.1;% Stopping criterion: relative decrease in the cost function smaller than 0.1
params.lambda         = 2^-22;

[ U, V ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% For example, we display dictionary elements 1,10 and 20

imagesc( reshape( V(:,1), Height, Width )' )

pause

imagesc( reshape( V(:,10), Height, Width )' )

pause

imagesc( reshape( V(:,20), Height, Width )' )

pause

%==========================================================================
%==========================================================================
% We keep the same setting as above, but we also show the effect of
% controlling the number of different nonzero patterns among the r
% dictionary elements.
% The parameter to set is params.m, such that mod(r,params.m)==0.

clear mode
mode.rectangle = true;
mode.pi4       = true;

spG = get_groups( mode, [ Height, Width ] );

%--------------------------------------------------------------------------
% We still consider the Lq-L2 quasi-norm regularization, q in (0,1)
clear params

q = 0.5;
params.normparam      = q;
params.r              = r;
params.it0            = 25;% Cost function displayed every 25 iterations
params.min_delta_cost = 0.1;% Stopping criterion: relative decrease in the cost function smaller than 0.1
params.lambda         = 2^-22;

%--------------------------------------------------------------------------
% There are 18 different nonzero patterns, each
% one represented by r/m=2 dictionary elements.
params.m = 18;

[ U, V ] = sspca( X, spG, params );

%--------------------------------------------------------------------------
% For example, we display dictionary elements 3,4 that share the same
% nonzero pattern

imagesc( reshape( V(:,3), Height, Width )' )

pause

imagesc( reshape( V(:,4), Height, Width )' )

pause

% We also display dictionary elements 9,10 that share the same
% nonzero pattern

imagesc( reshape( V(:,9), Height, Width )' )

pause

imagesc( reshape( V(:,10), Height, Width )' )

pause

close all

echo off
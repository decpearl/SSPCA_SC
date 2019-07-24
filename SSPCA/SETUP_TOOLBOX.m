path_to_src = '.';% *** TO BE FILLED ***

addpath(path_to_src);
addpath([path_to_src,'/Mex/']);
addpath([path_to_src,'/Misc/']);

% COMPILE MEX FILES

cd([path_to_src '/Mex'])

mex mexUpdateEta.cpp       -largeArrayDims
mex mexGetOmega.cpp        -largeArrayDims
mex mexUpdateU.cpp -l blas -largeArrayDims
mex mexUpdateV.cpp -l blas -largeArrayDims

cd ../

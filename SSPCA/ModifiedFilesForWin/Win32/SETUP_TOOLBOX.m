path_to_src = '.';% *** TO BE FILLED ***

addpath(path_to_src);
addpath([path_to_src,'/Mex/']);
addpath([path_to_src,'/Misc/']);

% COMPILE MEX FILES

cd([path_to_src '/Mex'])

mex mexUpdateEta.cpp       -largeArrayDims
mex mexGetOmega.cpp        -largeArrayDims

lapacklib = fullfile(matlabroot,'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
blaslib   = fullfile(matlabroot, 'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');

mex('-l', '-largeArrayDims', 'mexUpdateU.cpp', blaslib, lapacklib);
mex('-l', '-largeArrayDims', 'mexUpdateV.cpp', blaslib, lapacklib);

cd ../

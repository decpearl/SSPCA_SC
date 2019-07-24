
clc, clear all, close all

addpath('SSPCA');
addpath('SSPCA/Mex/');
addpath('SSPCA/Misc/'); 
addpath('ClusteringMeasure');


%read images
data_path = 'ORL\';

if ~isempty(rdir([data_path '**\*.pgm']))
    img_dir = rdir([data_path '**\*.pgm']);
    img_len = length(img_dir);
else
    return;
end
  

folder_string = data_path;
f = 0; n = 1;
%X = []; 
Y = []; s = []; featureVector = [];


for i=1:img_len
    disp(['file no.' int2str(i)])
    [pathstr, file_name, ext] = fileparts(img_dir(i).name);
    if strcmp(folder_string,pathstr) == 0,
        f = f+1;
        n = 1;
    else
        n = n+1;
    end
    folder_string = pathstr;
    % Read first frame.
    i_frame = im2double(imread(img_dir(i).name));
    [featureVector] = extractHOGFeatures(i_frame);
        
    Height = size(i_frame, 1);
    Width = size(i_frame, 2);
    Y = [Y; featureVector];
    s = [s f];
end

Y = im2double(Y);
n = size(unique(s),2);


e = size(Y,2);
a = round(e/10);

%------------------ SSPCA
 G            = zeros(e,10);  % training number = 20
    G(1:a,1)    = 1;
    G((a*1)+1:a*2,2)  = 1;
    G((a*2)+1:a*3,3) = 1;
    G((a*3)+1:a*4,4) = 1;
    G((a*4)+1:a*5,5) = 1;
    G((a*5)+1:a*6,6) = 1;
    G((a*6)+1:a*7,7) = 1;
    G((a*7)+1:a*8,8) = 1;
    G((a*8)+1:a*9,9) = 1;
    G((a*9)+1:e,10) = 1; 
spG          = sparse(G);

%--------------------------------------------------------------------------

acc = zeros(1,6);
k = 1;

    
spG          = sparse(G);
      
    clear params
    


q = 0.5;
params.normparam      = q;
params.r              = 50;%r=350 ; %r=350;  YaleB=350, Yale=550, ORL=450
params.it0            = 10;% Cost function displayed every 25 iterations
params.min_delta_cost = 0.0001;% Stopping criterion: relative decrease in the cost function smaller than 0.1
params.lambda         = 2^-22;
params.max_it         = 500;
params.intersection   = true;

[ U, V ] = sspca(Y, spG, params );

CMat = abs(U*U');

imagesc(CMat)


CKSym = CMat;

[CMatC,sc,OutlierIndx,Fail] = OutlierDetection(CMat,s);

gt = sc';


Grp = SpectralClustering(CKSym,n);
grps = bestMap(sc,Grp);


[A nmi avgent] = compute_nmi(gt,grps);
acc = Accuracy(grps,double(gt));
[f,p,r] = compute_f(gt,grps);
[ar,RI,MI,HI]=RandIndex(gt,grps); 
   



%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param = 1; % photon per pixel values to be simulated
ref_type = 'b'; %'b'=block, 's'=slit, 'p'=pinhole, 'a'=arbitrary
ref = []; % specify reference if ref_type='a'; otherwise can leave empty
%% Generate problem input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
save('data_in.mat')
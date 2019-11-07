% Generate the input file data_in.mat, which contains the variables:
% img - the input model ("image"))
% ref_type - the reference type ('b'-block, 'p'-pinhole, 's'-slit,
%                                                       'a'-arbitrary)
% ref - the reference values (can be empty array [] if ref_type is 'b','p',or 's')
% m - the length and width of the output data (assuming square output for now)
% photon_param - the specified number of photons per pixel
% (all other variables used in the process are cleared)
%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param = 100; % photon per pixel values to be simulated
ref_type = 'b'; %'b'=block, 's'=slit, 'p'=pinhole, 'a'=arbitrary
ref = []; % specify reference if ref_type='a'; otherwise can leave empty
%% Generate problem input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
%%%%
clear n; clear namestr; clear stanstr; clear X; clear X_0;
%%%%
save('data_in.mat')
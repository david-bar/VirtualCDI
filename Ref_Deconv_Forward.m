% Forward model, which takes in data_in.mat, and outputs
% the file data_out.mat, which contains the variables:
% y - the simulated data (i.e. the squared Fourier transform magnitude
%     measurements with Poisson shot noise
% img_rec - the recovered image
% err - the relative recovery error between img and img_rec
% (all other variables used in the process are cleared)
load('data_in.mat')
rng(1);
n=size(img,1); %assuming square input for now
%% DFT matrices and partial DFT matrices
F1=dftmtx(m);
F1c=F1(:,[end-(n-1)+1:end,1:n]);
F1cc=F1(:,1:n);
F2=dftmtx(m);
F2c=F2(:,[end-(2*n-1)+1:end,1:2*n]);
F2cc=F2(:,1:n);

%% Reference choice

if ref_type=='b';
    ref=ones(n);
end

if ref_type=='s';
    ref=zeros(n); ref(:,end)=ones(n,1);
end

if ref_type=='p';
    ref=zeros(n); ref(n,n)=1;
end

if ref_type=='a'
    MR=ref2mtrx(ref);
end

%% Image and Reference Composite
x = [img, ref];

%% Squared Fourier Transform Magnitude 
xpad = zeros(m);
xpad(1:n, 1:2*n) = x;
f = fft2(xpad);
y_clean = abs(f).^2; 

n_photon = photon_param*m^2;
y = sum(y_clean(:))*(n_photon)^(-1) * poissrnd( (n_photon/sum(y_clean(:)))*y_clean );
%% Get autocorrelation
rfull=(1/(m^2))*(F1c')*y*transpose((F2c'));
r = rfull(1:n,1:n);
%% Run Referenced Deconvolution Algorithm
z=img_recov(r, n, ref, ref_type);
%% Error
err = norm(img(:)-z(:))/norm(img(:));
%% Save Recovered Image
img_rec=real(z);
save('data_out.mat')
%%
clear n; clear F1; clear F1c; clear F1cc; clear F2; clear F2c; clear F2cc;
if ref_type=='a'
    clear MR;
end
clear x; clear xpad; clear f; clear y_clean; clear n_photon; clear rfull; clear r;
clear z;

% Function which recovers a specimen image from a reference r of size n x
% n, and of type ref_type ('b'-block, 's'-slit, 'p'-pinhole, 'a'-arbitrary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=img_recov(r, n, ref, ref_type);

%% Run algorithm
if ref_type=='b'
        k=tril(ones(n));
        z = k \ r / k';
end

if ref_type=='s'
        k=tril(ones(n));
        z = k \ r / eye(n)';
end

if ref_type=='p'
            z=r;            
end
if ref_type=='a'
    L=sparse(ref2mtrx(ref));
    b=r; b=b(:);
    z=L\b;
end
z=reshape(z,n,n);

end
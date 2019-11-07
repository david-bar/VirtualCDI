load('data_in.mat')
% The important input variables in data_in.mat are:
% img - the input model ("image"))
% ref_type - the reference type ('b'-block, 'p'-pinhole, 's'-slit,
%                                                       'a'-arbitrary)
% ref - the reference values (can be empty array [] if ref_type is 'b','p',or 's')
% m - the length and width of the output data (assuming square output for now)
% photon_param - the specified number of photons per pixel

% The important output variables in data_out are:
% img_rec - the recovered image
% err - the relative squared recovery error between img and img_rec
% exp_err - the (analytically calculated) expected value for err
% S - the reference scaling factor corresponding to ref


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
    % Scaling factor
    S_dec = 1/(m^2) + 2*((n-1)/(m^2))*(1-cos(2*pi*[0:m-1]/m));
    S=S_dec'*S_dec;
end

if ref_type=='s';
    ref=zeros(n); ref(:,end)=ones(n,1);
    % Scaling factor
    S_dec = 1/(m^2) + 2*((n-1)/(m^2))*(1-cos(2*pi*[0:m-1]/m));
    S_flat = (n/(m^2))*ones(1,m);
    S=S_dec'*S_flat;
end

if ref_type=='p';
    ref=zeros(n); ref(n,n)=1;
    % Scaling factor
    S_flat = (n/(m^2))*ones(1,m);
    S=S_flat'*S_flat;
end

if ref_type=='a'
    MR=ref2mtrx(ref);
    TR=(1/(m^2))*inv(MR)*(kron(F2cc',F1cc'));
    S=reshape(diag(TR'*TR),m,m);
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
%% Observed Error
err = (norm(img(:)-z(:))/norm(img(:)))^2;
%% Expected Error
c=n_photon/sum(y_clean(:));
exp_err=(1/c)*sum(S(:).*y_clean(:))/(norm(img(:))^2);
%% Save Recovered Image
img_rec=real(z);
save('data_out.mat')
%%
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
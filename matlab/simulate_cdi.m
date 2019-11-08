function simulate_cdi(model_fname, output_fname, opts)

if nargin == 0 test_simulate_cdi(); return; end;

if nargin < 3 opts = struct(); end;

dopts = struct();
dopts.m = 1024; % detector dimension 
dopts.photon_param = 100; % photon per pixel values to be simulated
dopts.ref_type = 'b'; %'b'=block, 's'=slit, 'p'=pinhole, 'a'=arbitrary
dopts.ref = []; % specify reference if ref_type='a'; otherwise can leave empty

doptnames = fieldnames(dopts);
for i = 1:length(doptnames)
    if ~isfield(opts, doptnames{i})
        opts = setfield(opts, doptnames{i}, getfield(dopts, doptnames{i}));
    end;
end;

model = load(model_fname, 'opts', 'img');
img = model.img;
ref_type = opts.ref_type;

rng(1);
n=size(img,1); %assuming square input for now
m=opts.m;
photon_param=opts.photon_param;
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
    ref = opts.ref
    % MR=ref2mtrx(ref);
end

%% Image and Reference Composite
x = [img, ref];

%% Squared Fourier Transform Magnitude 
xpad = zeros(m);
xpad(1:n, 1:2*n) = x;
f = fft2(xpad);
y_clean = abs(f).^2; 

n_photon = photon_param*m^2;
data = sum(y_clean(:))*(n_photon)^(-1) * poissrnd( (n_photon/sum(y_clean(:)))*y_clean );

% %% Get autocorrelation
% rfull=(1/(m^2))*(F1c')*y*transpose((F2c'));
% r = rfull(1:n,1:n);
% %% Run Referenced Deconvolution Algorithm
% z=img_recov(r, n, ref, ref_type);
% %% Error
% err = norm(img(:)-z(:))/norm(img(:));
% %% Save Recovered Image
% img_rec=real(z);
% %%
% clear n; clear F1; clear F1c; clear F1cc; clear F2; clear F2c; clear F2cc;
% if ref_type=='a'
%     clear MR;
% end
% clear x; clear xpad; clear f; clear y_clean; 
% clear n_photon; clear rfull; clear r;
% clear z;


save(output_fname, 'opts', 'model', 'data');


function test_simulate_cdi

mfiledir = fileparts(mfilename('fullpath'));
specimen_fname = [mfiledir filesep '..' filesep 'Models' filesep 'plos10.png'];
create_model(specimen_fname, 'test_model.mat')

out_fname = 'test_simulated.mat';
disp(['Writing to ' out_fname]);
simulate_cdi('test_model.mat', out_fname)

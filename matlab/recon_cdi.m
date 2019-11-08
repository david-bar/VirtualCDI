function recon_cdi(simulated_fname, recon_out_fname, opts)

if nargin == 0 test_recon_cdi(); return; end;

if nargin < 3 opts = struct(); end;

dopts = struct();
dopts.n = 64; % image and reference dimension
dopts.ref_type = 'b'; % image and reference dimension
dopts.ref = [];

doptnames = fieldnames(dopts);
for i = 1:length(doptnames)
    if ~isfield(opts, doptnames{i})
        opts = setfield(opts, doptnames{i}, getfield(dopts, doptnames{i}));
    end;
end;

simulated = load(simulated_fname, 'opts', 'data', 'model');

data = simulated.data;
model = simulated.model;
m = simulated.opts.m;
n = opts.n;
ref_type = opts.ref_type;

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
end

%% DFT matrices and partial DFT matrices
F1=dftmtx(m);
F1c=F1(:,[end-(n-1)+1:end,1:n]);
F1cc=F1(:,1:n);
F2=dftmtx(m);
F2c=F2(:,[end-(2*n-1)+1:end,1:2*n]);
F2cc=F2(:,1:n);

%% Get autocorrelation
rfull=(1/(m^2))*(F1c')*data*transpose((F2c'));
r = rfull(1:n,1:n);
%% Run Referenced Deconvolution Algorithm
z=img_recov(r, n, ref, ref_type);
% %% Error
% err = norm(img(:)-z(:))/norm(img(:));
%% Save Recovered Image
img_recon=real(z);

save(recon_out_fname, 'opts', 'data', 'model', 'img_recon');

% Function which recovers a specimen image from a reference r of size n x
% n, and of type ref_type ('b'-block, 's'-slit, 'p'-pinhole, 'a'-arbitrary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=img_recov(r, n, ref, ref_type)

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

function test_recon_cdi

mfiledir = fileparts(mfilename('fullpath'));
specimen_fname = [mfiledir filesep '..' filesep 'Models' filesep 'plos10.png'];
create_model(specimen_fname, 'test_model.mat');

simulated_fname = 'test_simulated.mat';
disp(['Writing to ' simulated_fname]);
simulate_cdi('test_model.mat', simulated_fname);

recon_fname = 'test_recon.mat';
model = load('test_model.mat');
simulated = load(simulated_fname);
n = model.opts.n;
ref_type = simulated.opts.ref_type;
disp(['Writing to ' recon_fname]);
recon_cdi(simulated_fname, recon_fname, struct('n', n, 'ref_type', ref_type));


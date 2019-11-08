function create_model(specimen_fname, model_output_fname, opts)

if nargin == 0 test_create_model(); return; end;

if nargin < 3 opts = struct(); end;

dopts = struct();
dopts.n = 64; % image and reference dimension

doptnames = fieldnames(dopts);
for i = 1:length(doptnames)
    if ~isfield(opts, doptnames{i})
        opts = setfield(opts, doptnames{i}, getfield(dopts, doptnames{i}));
    end;
end;

X = mat2gray(imread(specimen_fname));
X_0 = rgb2gray(X);
img = imresize(X_0,[opts.n,opts.n]);

save(model_output_fname, 'opts', 'img');


function test_create_model

mfiledir = fileparts(mfilename('fullpath'));
specimen_fname = [mfiledir filesep '..' filesep 'Models' filesep 'plos10.png'];

out_fname = 'test.mat';
disp(['Writing to ' out_fname]);
create_model(specimen_fname, out_fname);
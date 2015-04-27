% Constructs gradient, accounting for different temporal resolutions

% lGrad - each col is gradient wrt linear filter at each time t
% psGrad - postspike filter
% k1 - 1st kernel (uses spike data)
% k2 - 2nd kernel (uses spike data) - not used in lite.

function g = multisample_gradHess_structs_lite(basepars,n,k1,~,trainpars)

% Get dimensions
stimT = size(trainpars.lgrad,2);
%spikeT = size(trainpars.psbasisGrad,2);
spikeT = size(trainpars.psbasisGrad{1},2);
if (mod(spikeT,stimT) > 0)
   error('ERROR: stim resolution is not integer multiple of spike resolution!stimT=%d,spikeT=%d\n',stimT,spikeT);
end
f = spikeT/stimT;

if size(k1,1) ~= spikeT 
   error('ERROR: kernels are not at the spike resolution!\n');
end

%nLpars = size(trainpars.lgrad{1},1);
nLpars = size(trainpars.lgrad,1);
nPSpars = size(trainpars.psbasisGrad{1},1);

npars = 1 + nLpars + nPSpars;

% Calculate "block-summed" kernel to use for coarser time resolution
k1blocked = sum(reshape(k1(:,n),f,stimT),1)'; % stimT x 1

% Init output variables
g = zeros(npars,1);

offset = 0;
b_idx = offset+1;

offset = offset + 1;
L_idx = offset+1:offset+nLpars;

offset = offset + nLpars;
PS_idx = offset+1:offset+nPSpars;

% Get the relevant indices used for training: edoi: CHECK
relidx_c = get_train_idx(basepars);
relidx_f = get_train_idx(basepars,'spike');

% Copy relevant matrices
%lgradn = trainpars.lgrad{n}(:,relidx_c);
lgradn = trainpars.lgrad(:,relidx_c);
%psgradn = trainpars.psbasisGrad{trainpars.baseneuron_idx(n)}(:,relidx_f);
psgradn = trainpars.psbasisGrad{1}(:,relidx_f);

k1 = k1(relidx_f,n);
k1b = k1blocked(relidx_c);

% Calculate the gradient - essentiallly weight the gradients by the kernel k1
% but have to account for different resolutions

%-- Base rate
g(b_idx) = sum(k1b);
%-- Linear filters (or stimulus-dependent term)
g(L_idx) = lgradn*k1b;
%-- Postspike filters
g(PS_idx) = psgradn*k1;


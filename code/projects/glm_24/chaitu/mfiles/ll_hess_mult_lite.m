function W = ll_hess_mult_lite(Hinfo,Y,basepars,~,trainpars)

% Function that indirectly compute H*Y where H is the hessian of log likelihood.

% Hinfo.kernel2 is the kernel (Tx1) vector (weighting vector)

% The trick is that H = A*K*A' where K is diagonal matrix
% Computing H*Y conventionally takes:
%    O(nT) for computing Hinfo*K;
%    O(n^2*T) for computing (Hinfo*K) * Hinfo'
%    O(n^2*b) for computing H*Y
% so the dominant term is O(n^2*T)

% Indirectly it takes:
%   O(n*T) for computing Hinfo*K
%   O(n*T*b) for computing Hinfo'*Y
%   O(n*T*b) for computing (Hinfo*K) * (Hinfo'*Y)
% so the dominant term is O(n*T*b) which is much less if n is large

% It will be assumed that ther diagonal entries of the K matrix will be stored in Hinfo.kernel2

%-- Determine relevant idx
relidx_c = get_train_idx(basepars);
relidx_f = get_train_idx(basepars,'spike');

%nLpars = size(trainpars.lgrad{1},1);
nLpars = size(trainpars.lgrad,1);
npars_perneuron = get_npars(basepars,1);

idx = 1:npars_perneuron;

%keyboard

coarsemat = [ones(1,length(relidx_c)); trainpars.lgrad(:,relidx_c)];
finemat   = trainpars.psbasisGrad{1}(:,relidx_f);

%keyboard
W = -multres_product(coarsemat, finemat, Hinfo.kernel2(relidx_f), Y(idx(1:1+nLpars),:), Y(idx(1+nLpars+1:end),:));

if (isfield(basepars,'frozen_idx'))
   W(basepars.frozen_idx,:) = 0;
end

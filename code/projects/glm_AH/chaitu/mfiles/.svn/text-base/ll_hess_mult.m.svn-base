function W = ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars)

% Function that indirectly compute H*Y where H is the hessian of log
% likelihood

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


% It will be assumed that ther diagonal entries of the K matrix will be
% stored in Hinfo.kernel2

%fprintf('I am being called\n');

 
% Determine relevant idx
relidx_c = get_train_idx(basepars);
relidx_f = get_train_idx(basepars,'spike');
%eval('whos  relidx_*')
%keyboard

1;
b = size(Y,2);
Neff = size(trainpars.D,2);
nLpars = size(trainpars.lgrad{1},1);
if (isfield(basepars,'XsqK') && basepars.XsqK)
    %trainpars.lsqgrad{1} = stimpars.dt .*(1/stimpars.dt.*trainpars.lgrad{1}).^2;
    trainpars.lsqgrad{1} = stimpars.dt .* sqrt(((1/stimpars.dt.*trainpars.lgrad{1})-0.5).^2);    
    nLpars = nLpars + size(trainpars.lsqgrad{1},1);
end
    
if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
      trainpars.extgrad{1} = interpmtx(basepars.interpfilt,basepars.ext_timepts,basepars.maxt)';%sincinterpmtx(makeaxis(stimpars.dt,basepars.maxt),basepars.ext_timepts)';
      nLpars = nLpars + size(trainpars.extgrad{1},1);
end

npars_perneuron = get_npars(basepars,Neff);
nneurons = size(Hinfo.kernel2,2);

W = zeros(npars_perneuron*nneurons,b);

for j=1:nneurons
    
    
    1;
    
    j_idx = (j-1)*npars_perneuron+1:j*npars_perneuron;
    
    
%    if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
%        coarsemat = trainpars.extgrad{j}(:,relidx_c);
%        coarsemat = [coarsemat;ones(1,length(relidx_c)); trainpars.lgrad{j}(:,relidx_c)];
%    else
%    coarsemat = [];
%    coarsemat = [coarsemat;ones(1,length(relidx_c)); trainpars.lgrad{j}(:,relidx_c)];
     coarsemat = [ones(1,length(relidx_c)); trainpars.lgrad{j}(:,relidx_c)];
     % this is a large matrix (edoi)
%end
    

    
%    if (isfield(basepars,'XsqK') && basepars.XsqK && isfield(trainpars,'lsqgrad') && size(trainpars.lsqgrad{j},1)>0)
%        coarsemat = [coarsemat; trainpars.lsqgrad{j}(:,relidx_c)];
%    end
    
    if (size(trainpars.D,2) > 1)
        finemat = [trainpars.psbasisGrad{j}(:,relidx_f); vertcat(trainpars.cpbasisGrad{[1:j-1 j+1:Neff]}(:,relidx_f))];
    else
        finemat = [trainpars.psbasisGrad{j}(:,relidx_f)];
    end
    
    W(j_idx,:) = -multres_product(coarsemat,finemat,Hinfo.kernel2(relidx_f,j),Y(j_idx(1:1+nLpars),:),Y(j_idx(1+nLpars+1:end),:));
    if (isfield(basepars,'frozen_idx'))
        W(basepars.frozen_idx,:) = 0;
    end
    
    % Determine if there is any regularization
    regflag = isfield(basepars,'stimreg') && strcmp(basepars.stimreg,'L2');
    if (regflag)
        stimidx_j = get_pars_idx(basepars,j,Neff,'k');
        W(stimidx_j,:) = W(stimidx_j,:) + basepars.lambda_reg .* Y(stimidx_j,:);
    end
    
end
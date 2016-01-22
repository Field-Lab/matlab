% Constructs gradient and hessian, accounting for different temporal
% resolutions

% lGrad - each col is gradient wrt linear filter at each time t
% psGrad - "" "" postspike filter ""
% cpGrad - "" "" coupling filter "" 
% k1 - 1st kernel (uses spike data)
% k2 - 2nd kernel (uses spike data)

function [g H] = multisample_gradHess_structs_fixrk2(basepars,n,k1,k2,trainpars,cpGrad,frozen_idx)

cpFlag = (nargin > 4 && size(trainpars.D,2) > 1);

% Get dimensions

%stimT = size(trainpars.lgrad{1},2);
stimT = size(trainpars.lgrad{1}{1},1); % edoi
spikeT = size(trainpars.psbasisGrad{1},2);

if (cpFlag && size(cpGrad,2) ~= spikeT)
    error('ERROR: dimension mismatch between coupling and postspike gradients!\n');
end

if (mod(spikeT,stimT) > 0)
    error('ERROR: stim resolution is not integer multiple of spike resolution!stimT=%d,spikeT=%d\n',...
       stimT,spikeT);
end

f = spikeT/stimT;

if (size(k1,1) ~= spikeT || (nargout > 1 && size(k2,1) ~= spikeT))
    error('ERROR: kernels are not at the spike resolution!\n');
end
1;

extflag = (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts));
if (extflag)
    nExtpars = length(basepars.ext_timepts);
else
    nExtpars = 0;
end

nLpars = 2; %size(trainpars.lgrad{1},1); edoi
if (isfield(basepars,'XsqK') && basepars.XsqK)
    nLsqpars = nLpars;
else
    nLsqpars = 0;
end

nPSpars = size(trainpars.psbasisGrad{1},1);
if (cpFlag)
    nCPpars = size(cpGrad,1);
else
    nCPpars = 0;
end

npars = nExtpars + 1 + nLpars + nLsqpars + nPSpars + nCPpars;

% Calculate "block-summed" kernel to use for coarser time resolution
k1blocked = sum(reshape(k1(:,n),f,stimT),1)'; % stimT x 1

% Init output variables

g = zeros(npars,1);

offset = 0;
% Indices for each group of parameters
if (extflag)
    ext_idx = 1:nExtpars;
    offset = offset + nExtpars;
else
    ext_idx = [];
end

b_idx = [offset+1];
offset = offset + 1;
L_idx = offset+1:offset+nLpars;
offset = offset + nLpars;

if (nLsqpars>0)
    Lsq_idx = offset+1:offset+nLsqpars;
    offset = offset + nLsqpars;
else
    Lsq_idx = [];
end

PS_idx = offset+1:offset+nPSpars;
offset = offset+nPSpars;

if cpFlag
    CP_idx = offset+1:npars;
else
    CP_idx = [];
end

% Get the relevant indices used for training
relidx_c = get_train_idx(basepars);
relidx_f = get_train_idx(basepars,'spike');


% Copy relevant matrices
if (nExtpars > 0)
    egradn = trainpars.extgrad{n}(:,relidx_c);
end

%--
fprintf('multisample_gradHess_structs_fixrk2.m \n')
keyboard
% edoi
%lgradn = trainpars.lgrad{n}(:,relidx_c);
lgradn = [trainpars.lgrad{1}{1},trainpars.lgrad{1}{2}]';
%--

if (nLsqpars > 0)
   lsqgradn = trainpars.lsqgrad{n}(:,relidx_c);
end
psgradn = trainpars.psbasisGrad{trainpars.baseneuron_idx(n)}(:,relidx_f);
if (cpFlag)
   cpgradn = cpGrad(:,relidx_f);
end

k1 = k1(relidx_f,n);
k1b = k1blocked(relidx_c);
if (~isempty(k2))
    k2 = k2(relidx_f,n);
end

% Calculate the gradient - essentially weight the gradients by the kernel k1
% but have to account for different resolutions

% Extrinsic parameters
if (nExtpars > 0)
    g(ext_idx) = egradn*k1b;
end

% Base rate
g(b_idx) = sum(k1b); % edoi: check - nan?

% Linear filters
g(L_idx) = lgradn*k1b;

if (nLsqpars > 0)
   g(Lsq_idx) = lsqgradn*k1b;
end

% Postspike filters
g(PS_idx) = psgradn*k1;

% Coupling filters
if (cpFlag)
    g(CP_idx) = cpgradn*k1;
end

%% for Hessian

if (nargout > 1)
   
   H = zeros(npars,npars);
   k2b = sum(reshape(k2,f,length(relidx_c)),1)'; %% stimT x 1
   
   % Calculate upper triangle of most of the hessian (2nd term of product rule - the 1st term is
   % precalculated and stored (trainpars.hst))
   
   
   % first compute blocked version of PS gradient weighted by k2
   psk2weighted = mult_col(psgradn,k2);
   psk2blocked = reshape(sum(reshape(psk2weighted,nPSpars,f,length(relidx_c)),2),nPSpars,...
      length(relidx_c)); % nPSpars x stimT
   
   if (cpFlag)
      cpk2weighted = mult_col(cpgradn,k2);
      %cpk2weighted = (cpGrad.*repmat(k2',nCPpars,1));
      cpk2blocked = reshape(sum(reshape(cpk2weighted,nCPpars,f,length(relidx_c)),2),nCPpars,...
         length(relidx_c)); % nCPpars x stimT
   end
   
   if (nExtpars > 0)
      
      % Extrinsic rate - Extrinsic rate
      ext2weighted = mult_col(egradn,k2b);
      H(ext_idx,ext_idx) = ext2weighted*egradn';
      
      % Extrinsic rate - Base rate
      H(ext_idx,b_idx) = sum(ext2weighted,2);
      
      % Extrinsic rate - linear filter
      H(ext_idx,L_idx) = ext2weighted*lgradn';
      
      % Extrinsic rate - square filter
      if (nLsqpars > 0)
         H(ext_idx,Lsq_idx) = ext2weighted*lsqgradn';
      end
      
      % Extrinsic rate - PS filters
      H(ext_idx,PS_idx) = egradn*psk2blocked';
      
      % Extrinsic rate - CP filters
      if (cpFlag)
         H(ext_idx,CP_idx) = egradn*cpk2blocked';
      end
      
   end
   
   % Base rate - Base rate
   H(b_idx,b_idx) = sum(k2b); % 1 x 1
   
   % Base rate - Linear filters
   H(b_idx,L_idx) = (lgradn*k2b)'; % 1 x nLpars
   
   % base rate - Linear square filters
   if (nLsqpars > 0)
      H(b_idx,Lsq_idx) = (lsqgradn  * k2b); % 1 x nLsqpars
   end
   
   % Base rate - PS filters
   H(b_idx,PS_idx) = (psgradn*k2)'; % 1 x nPSpars
   
   % Base rate - CP filters
   if (cpFlag)
      H(b_idx,CP_idx) = (cpgradn*k2)'; % 1 x nCPpars
   end
   
   % Linear filter - Linear filter (heavy if n is large!)
   link2weighted = mult_col(lgradn,k2b);
   
   H(L_idx,L_idx) = link2weighted*lgradn';
   %H(L_idx,L_idx) = (trainpars.lgrad.*repmat(k2blocked',nLpars,1))*trainpars.lgrad';
   % nLpars x nLpars
   
   % Linear filter - Square filter
   if (nLsqpars>0)
      H(L_idx,Lsq_idx) = link2weighted*(lsqgradn)'; % nLpars x nLsqpars
   end
   
   % Linear filter - PS filter
   %psk2weighted = (trainpars.psbasisGrad.*repmat(k2',nPSpars,1));
   
   H(L_idx,PS_idx) = lgradn*psk2blocked';% nLpars x nPSpars
   
   % Linear filter - CP filter
   if (cpFlag)
      H(L_idx,CP_idx) = lgradn*cpk2blocked';% nLpars x nCPpars
   end
   
   if (nLsqpars>0)
      
      link2weighted = mult_col(lsqgradn,k2b);
      
      % Square filter - Square filter
      H(Lsq_idx,Lsq_idx) = link2weighted*(lsqgradn)'; % nLsqpars x nLsqpars
      
      % Square filter - PS filters
      H(Lsq_idx,PS_idx) = (lsqgradn)*psk2blocked'; % nLsqpars x nPSpars
      
      % Square filter - CP filters
      if (cpFlag)
         H(Lsq_idx,CP_idx) = (lsqgradn) * cpk2blocked'; % nLsqpars x nCPpars
      end
   end
   
   % PS filter - PS filter (heavy)
   H(PS_idx,PS_idx) = psk2weighted*psgradn'; % nPSpars x nPSpars;
   
   % PS filter - CP filter (heavy)
   if (cpFlag)
      H(PS_idx,CP_idx) = psk2weighted*cpgradn'; % nPSpars x nCPpars
   end
   
   % CP filter - CP filter (heavy)
   if (cpFlag)
      H(CP_idx,CP_idx) = cpk2weighted*cpgradn'; % nCPpars x nCPpars
   end
   
   
   % Symmetrize H
   H = triu(H) + triu(H)' - diag(diag(H));
   
end
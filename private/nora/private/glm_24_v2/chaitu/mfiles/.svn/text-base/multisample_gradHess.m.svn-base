% Constructs gradient and hessian, accounting for different temporal
% resolutions

% lGrad - each col is gradient wrt linear filter at each time t
% psGrad - "" "" postspike filter ""
% cpGrad - "" "" coupling filter "" 
% k1 - 1st kernel (uses spike data)
% k2 - 2nd kernel (uses spike data)

function [g H] = multisample_gradHess(k1,k2,lGrad,lsqGrad,psGrad,cpGrad)

cpFlag = (nargin > 4);
lsqFlag = ~isempty(lsqGrad);

% Get dimensions

stimT = size(lGrad,2);
spikeT = size(psGrad,2);

if (cpFlag && size(cpGrad,2) ~= spikeT)
    fprintf('ERROR: dimension mismatch between coupling and postspike gradients!\n');
end

if (mod(spikeT,stimT) > 0)
    fprintf('ERROR: stim resolution is not integer multiple of spike resolution!stimT=%d,spikeT=%d\n',stimT,spikeT);
end

f = spikeT/stimT;

if (length(k1) ~= spikeT || (nargout > 1 && length(k2) ~= spikeT))
    fprintf('ERROR: kernels are not at the spike resolution!\n');
end
1;
nLpars = size(lGrad,1);
if (lsqFlag)
    nLsqpars = size(lsqGrad,1); 
else
    nLsqpars = 0;
end
nPSpars = size(psGrad,1);
if (cpFlag)
    nCPpars = size(cpGrad,1);
else
    nCPpars = 0;
end

npars = 1 + nLpars + nPSpars + nCPpars;
1;
% Calculate "block-summed" kernel to use for coarser time resolution
k1blocked = sum(reshape(k1,f,stimT),1)'; % stimT x 1


% Init output variables

g = zeros(npars,1);

% Indices for each group of parameters
b_idx = [1];
L_idx = 1+1:1+nLpars;
if (lsqFlag)
    Lsq_idx = 1+nLpars+1:1+nLpars+nLsqpars;
else
    Lsq_idx = [];
end
PS_idx = 1+nLpars+nLsqpars+1:1+nLpars+nLsqars+nPSpars;
if cpFlag
    CP_idx = 1+nLpars+nLsqpars+nPSpars+1:npars;
else
    CP_idx = [];
end

% Calculate the gradient - essentiallly weight the gradients by the kernel k1
% but have to account for different resolutions

% Base rate
g(b_idx) = sum(k1blocked);

% Linear filters
g(L_idx) = lGrad*k1blocked;

% Square filters
g(Lsq_idx) = lsqGrad*k1blocked;

% Postspike filters
g(PS_idx) = psGrad*k1;

% Coupling filters
if (cpFlag)
    g(CP_idx) = cpGrad*k1;  
end

if (nargout > 1)
    
    H = zeros(npars,npars);
    k2blocked = sum(reshape(k2,f,stimT),1)'; %% stimT x 1
    
    % Calculate upper triangle of most of the hessian (2nd term of product rule - the 1st term is
    % precalculated and stored (trainpars.hst))
    
    % Base rate - Base rate
    H(b_idx,b_idx) = sum(k2blocked); % 1 x 1
    
    % Base rate - Linear filters
    H(b_idx,L_idx) = (lGrad*k2blocked)'; % 1 x nLpars
    
    % Base rate - Square filters
    if (lsqFlag)
        H(b_idx,Lsq_idx) = (lsqGrad*k2blocked)'; % 1 x nLsqpars
    end
    
    % Base rate - PS filters
    H(b_idx,PS_idx) = (psGrad*k2)'; % 1 x nPSpars
    
    % Base rate - CP filters
    if (cpFlag)
        H(b_idx,CP_idx) = (cpGrad*k2)'; % 1 x nCPpars
    end
    
    % Linear filter - Linear filter (heavy if n is large!)
    1;
    link2weighted = mult_col(lGrad,k2blocked);
    
    
    H(L_idx,L_idx) = link2weighted*lGrad';
    %H(L_idx,L_idx) = (lGrad.*repmat(k2blocked',nLpars,1))*lGrad'; % nLpars x nLpars
    
    % Linear filter - Square filter
    if (lsqFlag)
        H(L_idx,Lsq_idx) = link2weighted*lsqGrad'; % nLpars x nLsqpars
    end
    
    % Linear filter - PS filter
    % first compute blocked version of PS gradient weighted by k2
    psk2weighted = mult_col(psGrad,k2);
    %psk2weighted = (psGrad.*repmat(k2',nPSpars,1));
    psk2blocked = reshape(sum(reshape(psk2weighted,nPSpars,f,stimT),2),nPSpars,stimT); % nPSpars x stimT
    H(L_idx,PS_idx) = lGrad*psk2blocked';% nLpars x nPSpars
    
    % Linear filter - CP filter
    if (cpFlag)
        cpk2weighted = mult_col(cpGrad,k2);
        %cpk2weighted = (cpGrad.*repmat(k2',nCPpars,1));
        cpk2blocked = reshape(sum(reshape(cpk2weighted,nCPpars,f,stimT),2),nCPpars,stimT); % nCPpars x stimT
        H(L_idx,CP_idx) = lGrad*cpk2blocked';% nLpars x nCPpars
    end
    
    if (lsqFlag)
        
        link2weighted = mult_col(lsqGrad,k2blocked);
        
        % Square filter - Square filter
        H(Lsq_idx,Lsq_idx) = link2weighted*lsqGrad';
        
        % Square filter - PS filter
        H(Lsq_idx,PS_idx) = lsqGrad*psk2blocked';
        
        % Square filter - CP filter
        H(Lsq_idx,CP_idx) = lsqGrad*cpk2blocked';
        
    end
    
    
    
    % PS filter - PS filter (heavy)
    H(PS_idx,PS_idx) = psk2weighted*psGrad'; % nPSpars x nPSpars;
    
    % PS filter - CP filter (heavy)
    if (cpFlag)
        H(PS_idx,CP_idx) = psk2weighted*cpGrad'; % nPSpars x nCPpars
    end
    
    % CP filter - CP filter (heavy)
    if (cpFlag)
        H(CP_idx,CP_idx) = cpk2weighted*cpGrad'; % nCPpars x nCPpars
    end
    
    
    % Symmetrize H
    H = triu(H) + triu(H)' - diag(diag(H));
    
end
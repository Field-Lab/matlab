function W = ll_hess_mult_precon(Hinfo,Y,basepars,stimpars,trainpars,Mprecon)

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
1;
b = size(Y,2);
Neff = size(trainpars.D,2);
nLpars = size(trainpars.lgrad{1},1);
if (Neff > 1)
    nOpars = size(trainpars.psbasisGrad{1},1) + (Neff-1)*size(trainpars.cpbasisGrad{1},1);
else
    nOpars = size(trainpars.psbasisGrad{1},1);
end
Tcoarse = size(trainpars.lgrad{1},2);
npars_perneuron = 1+ nLpars + nOpars;
nneurons = size(Hinfo.kernel2,2);

W = zeros(npars_perneuron*nneurons,b);

nsamples = 5000;

%[sorted sorted_idx] = sort(abs(Hinfo.kernel2),'descend');

% Apply the preconditioner to Y

MtY = Mprecon'*Y;


for j=1:nneurons
    
    
    1;
    
    j_idx = (j-1)*npars_perneuron+1:j*npars_perneuron;
    
    if 0
        sample_idx = ceil(rand(1,nsamples).*size(stimpars.x,2));
        
        coarsemat = [ones(1,nsamples); trainpars.lgrad{j}(:,sample_idx)]; % usually just take all samples in order
        
        finesample_idx = (sample_idx-1).*basepars.fac; % starting indices
        finesample_idx = repmat((1:basepars.fac)',1,nsamples) + repmat(finesample_idx,basepars.fac,1);
        
        finesample_idx = finesample_idx(:);
        
        1;
        if (size(trainpars.D,2) > 1)
            %finemat = [trainpars.psbasisGrad{j}(:,finesample_idx);
            %vertcat(trainpars.cpbasisGrad{j}([1:j-1
            %j+n:Neff,finesample_idx)]; FIX IF YOU WANT TO USE
        else
            finemat = [trainpars.psbasisGrad{j}(:,finesample_idx)];
        end
        fprintf('Sampling %d timepoints to compute hessian matrix-vector product.\n',nsamples);
        
        W(j_idx,:) =-multres_product(coarsemat,finemat,Hinfo.kernel2(finesample_idx),MtY(1:1+nLpars,:),MtY(1+nLpars+1:end,:)).*(Tcoarse/nsamples);
        
    else
        
        coarsemat = [ones(1,size(trainpars.lgrad{1},2)); trainpars.lgrad{j}];
        if (size(trainpars.D,2) > 1)
            finemat = [trainpars.psbasisGrad{j}; vertcat(trainpars.cpbasisGrad{[1:j-1 j+1:Neff]})];
        else
            finemat = [trainpars.psbasisGrad{j}];
        end
        1;
        W(j_idx,:) = -Mprecon*multres_product(coarsemat,finemat,Hinfo.kernel2(:,j),MtY(j_idx(1:1+nLpars),:),MtY(j_idx(1+nLpars+1:end),:));
    end
    
end


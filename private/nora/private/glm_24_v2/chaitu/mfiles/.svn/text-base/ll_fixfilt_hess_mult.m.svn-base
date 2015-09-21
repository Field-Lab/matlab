function W = ll_fixfilt_hess_mult(Hinfo,Y,basepars,stimpars,trainpars)

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

b = size(Y,2);
nLpars = 1;
nOpars = size(trainpars.psbasisGrad,1) + (basepars.Nneurons - 1)*size(trainpars.cpbasisGrad,1);
Tcoarse = size(stimpars.x,2);
npars_perneuron = 1+ nLpars + nOpars;
nneurons = basepars.Nneurons;

W = zeros(npars_perneuron*nneurons,b);

for j=1:nneurons
    
    j_idx = (j-1)*npars_perneuron+1:j*npars_perneuron;
    
    coarsemat = [ones(1,Tcoarse); Hinfo.kx(:,j)'];
    %finemat = [trainpars.psbasisGrad(:,:,j); trainpars.cpbasisGrad(:,:,j)];  

    
    W(j_idx,:) = multres_product(coarsemat,trainpars.psbasisGrad(:,:,j),Hinfo.kernel2(:,j),Y(j_idx(1):j_idx(1+nLpars),:),Y(j_idx(1+nLpars+1:end),:));
    
end

W = -W;

if 0 % this code is not in use

    if (nargout > 2)
        
        
        
        
       idx = (j-1)*npartpars_perneuron:j*nartpars_perneuron;
       
       % Base rate - base rate
       Hinfo(idx(1),idx(1)) = sum(k2blocked(:,j));
       
       % Base rate - stim filt mag
       Hinfo(idx(1),idx(2)) = dot(kx(:,j),k2blocked);
       
       % Base rate - ps filters
       Hinfo(idx(1),idx(3:2+basepars.nofilters_postspike)) = basepars.psbasisGrad(:,:,j)*kernel2;
       
       % Base rate - cp filters
       if (basepars.Nneurons > 1)
           Hinfo(idx(1),idx(2+basepars.nofilters_postspike+1:2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling)) = cpGrad*kernel2;        
       end
       
       % stim filt mag - stim filt mag 
       Hinfo(idx(2),idx(2)) = dot(kx(:,j).*k2blocked,kx(:,j));
       
       % stim filt mag - ps filters
       psweighted_j = mult_col(basepars.psgradBasis(:,:,j),kernel2);
       Hinfo(idx(2),idx(3:2+basepars.nofilters_postspike)) = kx(:,j)' * block_spikes(psweighted_j,basepars.fac)';
       
       % stim filt mag - cp filters
       if (basepars.Nneurons > 1)
           cpweighted_j = mult_col(cpGrad,kernel2);
           Hinfo(idx(2),idx(2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling)) = kx(:,j)' * block_spikes(cpweighted_j,basepars.fac)';
       end
       
       % ps filters - ps filters
       Hinfo(idx(3:2+basepars.nofilters_postspike),idx(3:2+basepars.nofilters_postspike)) = psweighted_j*basepars.psgradBasis(:,:,j)';
       
       % ps filters - cp filters
       if (basepars.Nneurons > 1)
           Hinfo(idx(3:2+basepars.nofilters_postspike),idx(2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling)) = psweighted_j*cpGrad';
       end
       
       % cp filters - cp filters
       if (basepars.Nneurons > 1)
           Hinfo(idx(2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling),idx(2+basepars.nofilters_postspike+(basepars.Nneurons-1)*basepars.nofilters_coupling)) = cpweighted_j*cpGrad';
       end
        
    end

end
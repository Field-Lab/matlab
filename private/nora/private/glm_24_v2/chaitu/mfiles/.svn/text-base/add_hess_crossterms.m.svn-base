function ll_hess = add_hess_crossterms(ll_hess0,offset,basepars,stimpars,kernel,n)

1;

[nspace ntime] = get_nspacetime(basepars);

% Relevant idx for training
relidx = get_train_idx(basepars);

ll_hess = ll_hess0;
T = length(kernel);
% The only off-diagonal contributions to the hessian come from the
% linear filters

% Compute the missing term in the hessian wrt spatial and temporal filters (nspace x basepars.nofilters_k matrix)
[xidx1 yidx1 xidx2 yidx2] = get_sep_filt_idces(offset,basepars);

h = zeros(nspace,ntime);
switch(basepars.filtermode)
    case 'sep_basis'
        kernelblocked = sum(reshape(kernel(:,n),basepars.fac,T/basepars.fac),1)';
        % the hessian w.r.t. ith spatial and jth temporal basis fn
        % is the stimulus projected by the ith spatial filter and
        % convolved with the jth temporal filter.
        xproj = basepars.kspace_basis' * stimpars.x(basepars.crop_idx(:,n),:); % nspace x stimT
        
        
        for i=1:nspace
            for j=1:ntime
                1;
                proj_ij = stimpars.dt .* fastconv(xproj(i,:),basepars.ktime_basis(:,j)',1,size(stimpars.x,2));
                h(i,j) = dot(proj_ij(relidx),kernelblocked(relidx));                
            end
        end
    case 'sep_raw'
        kernelblocked = sum(reshape(kernel(:,n),basepars.fac,T/basepars.fac),1)';
        for k=1:ntime
            
            v = [basepars.padval.*ones(basepars.n,k-1) stimpars.x(basepars.crop_idx(:,n),1:(size(stimpars.x,2)-(k-1)))];
            
            
            h(:,k) = stimpars.dt .* v(:,relidx)*kernelblocked(relidx);
        end
end

ll_hess(xidx1,yidx1) = ll_hess(xidx1,yidx1) + h;
ll_hess(yidx1,xidx1) = ll_hess(yidx1,xidx1) + h';
ll_hess(xidx2,yidx2) = ll_hess(xidx2,yidx2) + h;
ll_hess(yidx2,xidx2) = ll_hess(yidx2,xidx2) + h';
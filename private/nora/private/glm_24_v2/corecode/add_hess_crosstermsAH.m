function ll_hess = add_hess_crosstermsAH(ll_hess0,offset,Basepars,Stimpars,kernel_grad,n)
%%%  called from train_ll_grad5
 %%%ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,Stimpars,kernel_grad,n);
 
k_space = Basepars.k_spacepixels;
k_time  = Basepars.k_stimframes;

% Relevant idx for training
relidx = get_train_idxAH(Basepars);

ll_hess = ll_hess0;
microbins = length(kernel_grad);
% The only off-diagonal contributions to the hessian come from the
% linear filters

% Compute the missing term in the hessian wrt spatial and temporal filters (k_space x Basepars.nofilters_k matrix)
[spaceidx1 timeidx1 spaceidx2 timeidx2] = get_sep_filt_idcesAH(offset,Basepars);

h = zeros(k_space,k_time);
switch(Basepars.k_filtermode)
    case 'sep_basis'
        kernel_gradblocked = sum(reshape(kernel_grad(:,n),Basepars.spikebins_perstimframe,T/Basepars.spikebins_perstimframe),1)';
        % the hessian w.r.t. ith spatial and jth temporal basis fn
        % is the stimulus projected by the ith spatial filter and
        % convolved with the jth temporal filter.
        xproj = Basepars.kspace_basis' * Stimpars.movie_ROI(Basepars.crop_idx(:,n),:); % k_space x stimT
        
        
        for i=1:k_space
            for j=1:k_time
                1;
                proj_ij = Stimpars.dt .* fastconv(xproj(i,:),Basepars.ktime_basis(:,j)',1,size(Stimpars.movie_ROI,2));
                h(i,j) = dot(proj_ij(relidx),kernel_gradblocked(relidx));                
            end
        end
    case 'sep_raw'
        kernel_gradblocked = sum(reshape(kernel_grad(:,n),Basepars.spikebins_perstimframe,microbins/Basepars.spikebins_perstimframe),1)';
        for k=1:k_time
            
            v = [Basepars.padval.*ones(Basepars.n,k-1) Stimpars.movie_ROI(Basepars.crop_idx(:,n),1:(size(Stimpars.movie_ROI,2)-(k-1)))];
            
            
            h(:,k) = Stimpars.dt .* v(:,relidx)*kernel_gradblocked(relidx);
        end
        
    case 'rk2'
        kernel_gradblocked = sum(reshape(kernel_grad(:,n),Basepars.spikebins_perstimframe,microbins/Basepars.spikebins_perstimframe),1)'; 
        
        
        
       %  ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,Stimpar
  %  case 'rk2'
  %      kernel_gradblocked = sum(reshape(kernel_grad(:,n),Basepars.spikebins_perstimframe,T/Basepars.spikebins_perstimframe),1)';
        % the hessian w.r.t. ith spatial and jth temporal basis fn
        % is the stimulus projected by the ith spatial filter and
        % convolved with the jth temporal filter.
  %      xproj = Basepars.kspace_basis' * Stimpars.movie_ROI(Basepars.crop_idx(:,n),:); % k_space x stimT
        
        
   %     for i=1:k_space
    %        for j=1:k_time
      %          1;
     %           proj_ij = Stimpars.dt .* fastconv(xproj(i,:),Basepars.ktime_basis(:,j)',1,size(Stimpars.movie_ROI,2));
      %          h(i,j) = dot(proj_ij(relidx),kernel_gradblocked(relidx));                
     %       end
    %    end
end

ll_hess(spaceidx1,timeidx1) = ll_hess(spaceidx1,timeidx1) + h;
ll_hess(timeidx1,spaceidx1) = ll_hess(timeidx1,spaceidx1) + h';
ll_hess(spaceidx2,timeidx2) = ll_hess(spaceidx2,timeidx2) + h;
ll_hess(timeidx2,spaceidx2) = ll_hess(timeidx2,spaceidx2) + h';
function H = hessian_sep_normconstr(x,lambda,basepars,stimpars,trainpars)

if (~(strcmp(basepars.filtermode,'sep_raw') || strcmp(basepars.filtermode,'sep_basis') ))
    fprintf('ERROR: hessian_sep_normconstr called in nonsep mode.\n');
end

[f g H] = ll_func2(x,basepars,stimpars,trainpars);
%npars_perneuron = get_npars(basepars,size(trainpars.D,2));

for j=1:basepars.Nneurons
    
    offset = get_pars_idx(basepars,j,size(trainpars.D,2),'k');
    offset = offset(1)-1;
    [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(offset,basepars);
    
    H(stim_idx1,stim_idx1) = H(stim_idx1,stim_idx1) + lambda.ineqnonlin(2*j-1).*eye(length(stim_idx1));
    H(stim_idx2,stim_idx2) = H(stim_idx2,stim_idx2) + lambda.ineqnonlin(2*j).*eye(length(stim_idx2));
    
    if (isfield(basepars,'XsqK') && basepars.XsqK)
        offset = get_pars_idx(basepars,j,size(trainpars.D,2),'ksq');
        offset = offset(1)-1;
        [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(offset,basepars);
    
        H(stim_idx1,stim_idx1) = H(stim_idx1,stim_idx1) + lambda.ineqnonlin(2*basepars.Nneurons+2*j-1).*eye(length(stim_idx1));
        H(stim_idx2,stim_idx2) = H(stim_idx2,stim_idx2) + lambda.ineqnonlin(2*basepars.Nneurons+2*j).*eye(length(stim_idx2));
    end
    
end


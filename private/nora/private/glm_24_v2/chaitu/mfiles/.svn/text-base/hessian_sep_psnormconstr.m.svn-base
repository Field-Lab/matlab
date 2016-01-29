function H = hessian_sep_psnormconstr(x,lambda,basepars,stimpars,trainpars)

if (~(strcmp(basepars.filtermode,'sep_raw') || strcmp(basepars.filtermode,'sep_basis') ))
    fprintf('ERROR: hessian_sep_normconstr called in nonsep mode.\n');
end

[f g H] = ll_func2(x,basepars,stimpars,trainpars);
%npars_perneuron = get_npars(basepars,size(trainpars.D,2));

for j=1:basepars.Nneurons
    
    ps_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'ps');    
    H(ps_idx,ps_idx) = H(ps_idx,ps_idx) + lambda.ineqnonlin(j).*eye(length(ps_idx));    
    
end


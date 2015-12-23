function [c ceq GC GCeq] = norm_constr(x,basepars,stimpars,trainpars)

if (~isfield(basepars,'normconstr_p'))
    basepars.normconstr_p = 2; %p  in Lp norm constraint
end
p = basepars.normconstr_p;

ceq = 0;
GCeq = zeros(size(x));
%npars_perneuron = get_npars(basepars,size(trainpars.D,2));

sqflag = (isfield(basepars,'XsqK') && basepars.XsqK);

switch(basepars.filtermode)
    case 'nonsep'
        c = zeros((sqflag+1)*basepars.Nneurons,1);
        GC = sparse(zeros(size(x,1),(sqflag+1)*basepars.Nneurons));
    case {'sep_raw','sep_basis'}
        c = zeros((sqflag+1)*2*basepars.Nneurons,1);
        GC = sparse(zeros(size(x,1),(sqflag+1)*2*basepars.Nneurons));
end



switch(basepars.filtermode)
    case 'nonsep'
        for j=1:basepars.Nneurons
            stim_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'k');
            c(j) = 1/p * (sum(abs(x(stim_idx)).^p) - basepars.normconstr^p);
            GC(stim_idx,j) = abs(x(stim_idx)).^(p-1).*sign(x(stim_idx));
            fprintf('Lp norm of filter is %0.3f. Constraint is %0.3f.\n',sum(abs(x(stim_idx)).^p),basepars.normconstr^p);
        end
        
        if (sqflag)
            for j=1:basepars.Nneurons
                stim_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'ksq');
                c(basepars.Nneurons+j) = 1/p * (sum(abs(x(stim_idx)).^p) - basepars.normconstr^p);
                GC(stim_idx,basepars.Nneurons+j) = abs(x(stim_idx)).^(p-1).*sign(x(stim_idx));
                fprintf('Lp norm of filter is %0.3f. Constraint is %0.3f.\n',sum(abs(x(stim_idx)).^p),basepars.normconstr^p);
            end
        end
        
        
    case {'sep_raw','sep_basis'}
        for j=1:basepars.Nneurons % Constrain only the spatial filters
            stim_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'k');
            [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(stim_idx(1)-1,basepars);
            c(2*j-1) = 1/p * (sum(abs(x(stim_idx1)).^p) - basepars.normconstr^p);
            c(2*j) = 1/p * (sum(abs(x(stim_idx2)).^p) - basepars.normconstr^p);
            GC(stim_idx1,2*j-1) = abs(x(stim_idx1)).^(p-1) .* sign(x(stim_idx1));
            GC(stim_idx2,2*j) = abs(x(stim_idx2)).^(p-1) .* sign(x(stim_idx2));
        end
        1;
        
        if (sqflag)
            for j=1:basepars.Nneurons % Constrain only the spatial filters
                stim_idx = get_pars_idx(basepars,j,size(trainpars.D,2),'ksq');
                [stim_idx1 blah stim_idx2 blah] = get_sep_filt_idces(stim_idx(1)-1,basepars);
                c(2*basepars.Nneurons+2*j-1) = 1/p * (sum(abs(x(stim_idx1)).^p) - basepars.normconstr^p);
                c(2*basepars.Nneurons+2*j) = 1/p * (sum(abs(x(stim_idx2)).^p) - basepars.normconstr^p);
                GC(stim_idx1,2*basepars.Nneurons+2*j-1) = abs(x(stim_idx1)).^(p-1) .* sign(x(stim_idx1));
                GC(stim_idx2,2*basepars.Nneurons+2*j) = abs(x(stim_idx2)).^(p-1) .* sign(x(stim_idx2));
            end
            
            
        end

        
        
end


c'
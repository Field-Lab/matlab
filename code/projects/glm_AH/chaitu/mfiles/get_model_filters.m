% Extract linear, postspike, and coupling filters from a pars GLM struct
% for neuron j
function [K PS C KSQ] = get_model_filters(p,j,pars,neff)

K = zeros(pars.n,pars.Mk);

if (nargout > 1)
   PS = zeros(pars.Mhist,1);
end

if (nargout > 2)
   C = zeros(pars.Mcoup,neff-1);
end


%npars = get_npars(pars,neff);

offset = get_pars_idx(pars,j,neff,'k');
offset = offset(1)-1;

switch pars.filtermode
   
   case 'sep_basis'
      
      [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,pars);
      K = (pars.kspace_basis*p(s1_idx))*(pars.ktime_basis*p(t1_idx))' + (pars.kspace_basis*p(s2_idx))*(pars.ktime_basis*p(t2_idx))';
      
   case {'sep_raw','rk2'}
      
      [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,pars);
      K = p(s1_idx)*p(t1_idx)' + p(s2_idx)*p(t2_idx)';
      
   case {'nonsep','raw'}
      
      K = reshape(p(offset+1:offset+pars.n*pars.Mk),pars.n,pars.Mk);
      
   case 'fixfilt'
      % commented out, edoi, 2012-01-05
      %K = pars.K{j}*p(offset+1);
end

if (nargout > 1)
    PS = pars.postspike_basis*p(get_pars_idx(pars,1,neff,'ps'));
else
    PS = [];
end

if (nargout > 2 && neff > 1)
    
    C = pars.coupling_basis*reshape(p(get_pars_idx(pars,1,neff,'cp')),pars.nofilters_coupling,(neff-1));
else
    C =[];
end

if (nargout > 3)
    if (~isfield(pars,'XsqK') || ~pars.XsqK)
        KSQ = [];
        return;
    end
    idx = get_pars_idx(pars,j,neff,'ksq');
    offset = idx(1)-1;
    switch pars.filtermode
        
        case 'sep_basis'
            
            [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,pars);
            KSQ = (pars.kspace_basis*p(s1_idx))*(pars.ktime_basis*p(t1_idx))' + (pars.kspace_basis*p(s2_idx))*(pars.ktime_basis*p(t2_idx))';
            
        case 'sep_raw'
            
            [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,pars);
            KSQ = p(s1_idx)*p(t1_idx)' + p(s2_idx)*p(t2_idx)';
            
        case 'nonsep'
            
            KSQ = reshape(p(offset+1:offset+pars.n*pars.Mk),pars.n,pars.Mk);

        case 'fixfilt'
            KSQ = pars.K{j}*p(offset+1);
    end
end
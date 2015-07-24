% Extract linear, postspike, and coupling filters from a pars GLM struct
% for neuron j
function [K PS C KSQ] = get_model_filtersAH(p,j,Basepars,neff)

K = zeros(Basepars.k_spacepixels,Basepars.k_stimframes);

if (nargout > 1)
   PS = zeros(Basepars.ps_timebins,1);
end

if (nargout > 2)
   C = zeros(Basepars.Mcoup,neff-1);
end


%npars = get_npars(pars,neff);

offset = get_pars_idxAH(Basepars,j,neff,'k');
offset = offset(1)-1;

switch Basepars.k_filtermode
   
   case 'sep_basis'
      
      [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idcesAH(offset,Basepars);
      K = (Basepars.kspace_basis*p(s1_idx))*(Basepars.ktime_basis*p(t1_idx))' + (Basepars.kspace_basis*p(s2_idx))*(Basepars.ktime_basis*p(t2_idx))';
      
   case {'sep_raw','rk2'}
      
      [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idcesAH(offset,Basepars);
      K = p(s1_idx)*p(t1_idx)' + p(s2_idx)*p(t2_idx)';
      
   case {'nonsep','raw'}
      
      K = reshape(p(offset+1:offset+Basepars.k_spacepixels*Basepars.k_stimframes),Basepars.k_spacepixels,Basepars.k_stimframes);
      
   case 'fixfilt'
      % commented out, edoi, 2012-01-05
      %K = Basepars.K{j}*p(offset+1);
end

if (nargout > 1)
    PS = Basepars.ps_basis*p(Basepars.paramind.PS);
else
    PS = [];
end

if (nargout > 2 && neff > 1)
    
    C = Basepars.coupling_basis*reshape(p(get_pars_idxAH(Basepars,1,neff,'cp')),Basepars.cp_filternumber,(neff-1));
else
    C =[];
end

if (nargout > 3)
    if (~isfield(Basepars,'XsqK') || ~Basepars.XsqK)
        KSQ = [];
        return;
    end
    idx = get_pars_idxAH(Basepars,j,neff,'ksq');
    offset = idx(1)-1;
    switch Basepars.k_filtermode
        
        case 'sep_basis'
            
            [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idcesAH(offset,Basepars);
            KSQ = (Basepars.kspace_basis*p(s1_idx))*(Basepars.ktime_basis*p(t1_idx))' + (Basepars.kspace_basis*p(s2_idx))*(Basepars.ktime_basis*p(t2_idx))';
            
        case 'sep_raw'
            
            [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,Basepars);
            KSQ = p(s1_idx)*p(t1_idx)' + p(s2_idx)*p(t2_idx)';
            
        case 'nonsep'
            
            KSQ = reshape(p(offset+1:offset+Basepars.k_spacepixels*Basepars.k_stimframes),Basepars.k_spacepixels,Basepars.k_stimframes);

        case 'fixfilt'
            KSQ = Basepars.K{j}*p(offset+1);
    end
end
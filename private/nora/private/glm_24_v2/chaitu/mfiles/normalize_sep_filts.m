function pstar2 = normalize_sep_filts(pstar,basepars,Neff,lp)

if (nargin < 4)
    lp = 2;
end

pstar2 = pstar;

if (~(strcmp(basepars.filtermode,'sep_raw') || strcmp(basepars.filtermode,'sep_basis')))
    return;
end

for j=1:basepars.Nneurons
    
    offset = get_pars_idx(basepars,j,Neff,'k');
    offset = offset(1)-1;
    [s1 t1 s2 t2] = get_sep_filt_idces(offset,basepars);
    
    n = sum(abs(pstar(s1)).^lp)^(1/lp); %Lp norm of 1st spatial filter
    pstar2(s1) = pstar2(s1)./n; %normalize 
    pstar(t1) = pstar(t1).*n; % absorb in temporal filter
    
    n = sum(abs(pstar(s2)).^lp)^(1/lp); % Lp norm of 2nd spatial filter
    pstar2(s2) = pstar2(s2)./n; %normalize 
    pstar(t2) = pstar(t2).*n; % absorb in temporal filter
    
    if (isfield(basepars,'XsqK') && basepars.XsqK)
        offset = get_pars_idx(basepars,j,Neff,'ksq');
        offset = offset(1)-1;
        [s1 t1 s2 t2] = get_sep_filt_idces(offset,basepars);
        
        n = sum(abs(pstar(s1)).^lp)^(1/lp); %Lp norm of 1st spatial filter
        
        1;
        
        pstar2(s1) = pstar2(s1)./n; %normalize
        pstar(t1) = pstar(t1).*n; % absorb in temporal filter
        
        n = sum(abs(pstar(s2)).^lp)^(1/lp); % Lp norm of 2nd spatial filter
        pstar2(s2) = pstar2(s2)./n; %normalize
        pstar(t2) = pstar(t2).*n; % absorb in temporal filter
    end
    
    
    
end
    
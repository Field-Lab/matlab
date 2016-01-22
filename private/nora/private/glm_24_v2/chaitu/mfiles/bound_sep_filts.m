function pstar2 = bound_sep_filts(pstar,basepars,Neff)

pstar2 = pstar;

if (~(strcmp(basepars.filtermode,'sep_raw') || strcmp(basepars.filtermode,'sep_basis')))
    return;
end

for j=1:basepars.Nneurons
    
    offset = get_pars_idx(basepars,j,Neff,'k');
    offset = offset(1)-1;
    [s1 t1 s2 t2] = get_sep_filt_idces(offset,basepars);
    
    n = max(abs(pstar(s1)));
    pstar2(s1) = pstar2(s1)./n; %normalize 
    pstar(t1) = pstar(t1).*n; % absorb in temporal filter
    
    n = max(abs(pstar(s2)));
    pstar2(s2) = pstar2(s2)./n; %normalize 
    pstar(t2) = pstar(t2).*n; % absorb in temporal filter
    
end
    
function pstar2 = adjust_baserate_ext(pstar,basepars,stimpars,trainpars)

pstar2 = pstar;

if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
    bold = pstar(get_pars_idx(basepars,1,size(trainpars.D,2),'b'));
    ext = pstar(get_pars_idx(basepars,1,size(trainpars.D,2),'ext'));
    [t ext_x] = get_ext_signal(ext,basepars,stimpars.dt);
    fprintf('Changing base firing rate from %0.3f to %0.3f\n',bold,bold + mean(ext_x));
    pstar2(get_pars_idx(basepars,1,size(trainpars.D,2),'b')) = bold + mean(ext_x);
    pstar2(get_pars_idx(basepars,1,size(trainpars.D,2),'ext')) = ext - mean(ext);
end
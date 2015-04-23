function pstar2 = sep2nonsep(pstar1,basepars)

K = get_model_filters(pstar1,1,basepars,1);

pstar2 = [pstar1(get_pars_idx(basepars,1,1,'ext')); ...
          pstar1(get_pars_idx(basepars,1,1,'b'));...
          K(:);...
          pstar1(get_pars_idx(basepars,1,1,'ps'));];

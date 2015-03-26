function init_pars = postspike_coupling_filterparamsAH(init_pars);
% MY_FUNCTION     This finds all the params associated with the 2 filters.
%                 Whcih were not specified in   
%   Calls creat_histbasis  (should be renamed create_psbasis)
%
%                                           'bore' - activate remotely
%  NOT A VERY GNERAL FUNCTION
%
% Works as of 
%



% Done in an microbin that's been adjusted for with the ps_fratio factor
%ps_f ration mean spost spike filter to linear filter ratio!!!
%init_pars.ps_timebins = floor((init_pars.spikebins_perstimframe*init_pars.k_stimframes)*init_pars.ps_fratio); % on fine timescale ~ 6/factor m(s)^-1




init_pars.ps_timebins = floor( (((init_pars.ps_ms)/1000) /init_pars.tstim) * init_pars.spikebins_perstimframe);
% linear filter basis -- not used.
init_pars.k_timebasis = []; init_pars.k_spacebasis = [];

% ps basis parameters     AH  literally just parameters or basis
% construction
train_dt= init_pars.tstim / init_pars.spikebins_perstimframe;
init_pars.ps_beta  = (init_pars.ps_timebins-1)*train_dt; % ending point for ps filters

% History filters (postspike and coupling)
bstretch = init_pars.ps_bstretch;
alpha_ps = init_pars.ps_alpha;
beta_ps  = init_pars.ps_beta;


if (init_pars.ps_filternumber > 0)
   init_pars.ps_basis = create_histbasis(alpha_ps,beta_ps,bstretch,...
      init_pars.ps_filternumber,init_pars.ps_timebins,train_dt,'L2',init_pars.ps_spacing);
else
   init_pars.ps_basis = [];
end
%%%%% Coupling filters%%%%%%%%%%%%%%%%%%
init_pars.cp_basis    =  [];

if init_pars.Coupling
   init_pars.cp_timebins = floor( (((init_pars.cp_ms)/1000) /init_pars.tstim) * init_pars.spikebins_perstimframe); % on fine timescale ~ 6/factor m(s)^-1


    % linear filter basis -- not used.
    init_pars.k_timebasis = []; init_pars.k_spacebasis = [];

    % cp basis parameters     AH  literally just parameters or basis
    % construction
    train_dt= init_pars.tstim / init_pars.spikebins_perstimframe;
    init_pars.cp_beta  = (init_pars.cp_timebins-1)*train_dt; % ending point for cp filters

    % History filters (postspike and coupling)
    bstretch = init_pars.cp_bstretch;
    alpha_cp = init_pars.cp_alpha;
    beta_cp  = init_pars.cp_beta;


    if (init_pars.cp_filternumber > 0)
       init_pars.cp_basis = create_histbasis(alpha_cp,beta_cp,bstretch,...
          init_pars.cp_filternumber,init_pars.cp_timebins,train_dt,'L2',init_pars.cp_spacing);
    else
       init_pars.cp_basis = [];
    end
end

if init_pars.Coupling && init_pars.BiDirect_CP
    
    CP_part1old                 = init_pars.cp_basis;
    CP_part1                    = [zeros( size (CP_part1old) ) ; CP_part1old];
    CP_part2                    = flipud (CP_part1);
    init_pars.cp_basis          = [CP_part1 , CP_part2];
end
    
    
    
    


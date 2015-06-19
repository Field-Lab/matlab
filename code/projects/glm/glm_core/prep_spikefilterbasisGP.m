% AKHeitman 2014-05-03
% Parameter DEPENDENT! GLMParams.spikefilter
% Cleaner version of creating spike filter basis previously used
% Calls creat_histbasis
function basis_vectors = prep_spikefilterbasisGP(basis_params,bin_size)
% MY_FUNCTION     This finds all the params associated with the 2 filters.
%                 Whcih were not specified in   
%   Calls create_histbasis  (should be renamed create_psbasis)
%
%                                           'bore' - activate remotely
%  NOT A VERY GNERAL FUNCTION
%
% Works as of 
%



% Done in an microbin that's been adjusted for with the ps_fratio factor
%ps_f ration mean spost spike filter to linear filter ratio!!!
%init_pars.ps_timebins = floor((init_pars.spikebins_perstimframe*init_pars.k_stimframes)*init_pars.ps_fratio); % on fine timescale ~ 6/factor m(s)^-1

dt = bin_size;
 
basis_params.timebins  = floor( (basis_params.ms/1000) / dt );
basis_params.beta    = (basis_params.timebins-1)*dt; % ending point for ps filters
% History filters (postspike and coupling)
bstretch = basis_params.bstretch;
alpha_ps = basis_params.alpha;
beta_ps  = basis_params.beta;
if (basis_params.filternumber > 0)
   basis_vectors = create_histbasis(alpha_ps,beta_ps,bstretch,...
      basis_params.filternumber,basis_params.timebins,dt,'L2',basis_params.spacing);
else
   basis_params.basis = [];
end

    
    


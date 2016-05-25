function  [f grad Hess_Base log_cif]= glm_opt_NonConvex_function(linear_params,covariates,lcif_derivative,spikebins,bin_duration)
 %%% PURPOSE %%%
% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized, column vector
% Covariates: time dependent input with multiplies into the params
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 

% Hess_Base still lacks the original correction term

% AKHEITMAN 2014-12-04
% NonConvex alteration done 2015-01-17
% For nonconvex   COVARIATE ~= LCIF_Derivative 
% because Covariates are dependent on the params

%%
% INITIALIZE PARAMS
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;
lcif_derivate = lcif_derivative;

% CONDITIONAL INTENSITY FUNCTIONS
lcif = p' * COV;
cif  = exp(lcif);


% THE OBJECTIVE FUNCTION (MONOTONIC IN LOG-LIKELIHOOD)
f_eval = sum( lcif(spt) ) - dt * sum(cif);

% GRADIENT FORMULA WITH LCIF_DERIVATIVE (RATHER THAN COV, DUE TO NONCONVEX)
g_eval = sum(lcif_derivate(:,spt),2)  - dt * ( lcif_derivate * (cif') );


% HESSIAN FORMULA WITH LCIF_DERIVATIVE (RATHER THAN COV, DUE TO NONCONVEX)
hessbase = zeros(size(lcif_derivate));
for i_vec = 1:size(lcif_derivate,1)
    hessbase(i_vec,:) = sqrt(cif) .* lcif_derivate(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');

% ADD PENALTY
penalty_strength = 500;
%L2
% f_penalty = penalty_strength*(p'*p);
% del_penalty = 2*penalty_strength*p;
% saturating
f_penalty = penalty_strength * sum(atan(p).^2);
del_penalty = 2*penalty_strength*atan(p).*(p.^2 + 1).^(-1);

% SWITCH SIGNS (TURN MAX TO MIN SO WE CAN USE FMINUNC OF MATLAB)
f            = -f_eval+f_penalty;
grad         = -g_eval+del_penalty;
Hess_Base    = -H_eval;
log_cif      = lcif;
end

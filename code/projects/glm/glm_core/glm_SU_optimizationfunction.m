function  [f grad Hess log_cif]= glm_SU_optimizationfunction(linear_params,covariates,spikebins,bin_duration)
 %%% PURPOSE %%%
% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% This one is designed for each linear_param to be a subunit weight (eg
%%% go through a nonlinearity)

%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized, column vector
% Covariates: time dependent input with multiplies into the params
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 


% AKHEITMAN 2014-12-04
%%
% Initialize
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;

f_eval = 0;
% Find Conditional Intensity and its log
for i = 1:length(p)
    lcif = p * COV(i,:);
    cif  = exp(lcif);
    
    % Evaluate the objective function (monotonic in log-likelihood)
    f_eval = f_eval + sum( lcif(spt) ) - dt * sum(cif);
end

% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );

% Evaluate the hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');


% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end
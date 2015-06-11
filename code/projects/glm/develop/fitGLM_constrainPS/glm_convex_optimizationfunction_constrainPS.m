function  [f grad Hess log_cif ps_filter]= glm_convex_optimizationfunction_constrainPS(linear_params,covariates,spikebins,bin_duration, ps_paramind,ps_basis,ps_balance)
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


% AKHEITMAN 2014-12-04

ps_filter0 = ps_basis * linear_params(ps_paramind);
drive = geomean(exp(ps_filter0));
if drive > 1
    ps_balance_coeff = log(1 / drive);
else
    ps_balance_coeff = 0;
end
ps_filter = ps_filter0 + ps_balance_coeff;

%%
% Initialize
p_adjust   = [linear_params; ps_balance_coeff];
COV_adjust = [covariates ; ps_balance];

% Find Conditional Intensity and its log
lcif = p_adjust' * COV_adjust;
cif  = exp(lcif);


COV = [covariates];
dt  = bin_duration;
spt = spikebins;
% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);

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
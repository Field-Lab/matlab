function  [f grad Hess log_cif]= glm_convex_optimizationfunction_squared(linear_params,covariates,spikebins,bin_duration)
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
%%

eps = 10^-5; 

% Initialize
p = linear_params;
COV = covariates;
dt = bin_duration;
% spt = spikebins;


% Find Conditional Intensity and its log
gen_signal = p' * COV;

% cutoff times
above_cutoff = (gen_signal > 0);
spt = zeros(size(above_cutoff));
spt(spikebins) = 1;
spt = logical(spt);

% F
cif = gen_signal.^2;
cif(gen_signal < 0) = eps;

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( log(cif(spt)) ) - dt * sum(cif);

% Evaluate the gradient
% g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );
g_eval = 2*sum(COV(:,spt & above_cutoff)./(repmat(gen_signal(spt & above_cutoff), [51,1])),2)-dt*2*sum(COV(:, above_cutoff).*repmat(gen_signal(above_cutoff), [51,1]),2);

% Evaluate the hessian
H_eval = -2*(COV(:,spt & above_cutoff)./repmat(gen_signal(spt & above_cutoff), [51,1]))*(COV(:,spt & above_cutoff)./repmat(gen_signal(spt & above_cutoff), [51,1]))' - dt*2*COV(:, above_cutoff)*COV(:, above_cutoff)';


% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = gen_signal;
end
function  [f grad Hess log_cif] = glm_SU_optimizationfunction(SU_params, SU_covariates, pooling_weights, time_filter, spikebins, bin_duration, non_stim_lcif)
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


% NB 2015-07-06

%% Initialize
dt = bin_duration;
spt = spikebins;
n_bins = size(non_stim_lcif);
n_params = length(SU_params);
penalty_strength = 0;

%% Find Conditional Intensity and its log
pixels = SU_covariates;
for i = 1:n_params
    SU_covariates(i, :, :) = SU_covariates(i, :, :) * SU_params(i);
end
stim_lcif = pooling_weights'* exp(squeeze(sum(SU_covariates,1)));
stim_lcif = conv(stim_lcif, flip(time_filter), 'same');
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;


% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));

% Add penalty
f_penalty = penalty_strength*SU_params'*SU_params;

%% Evaluate the gradient
% gradient of the lcif
del_LL = zeros(n_params, 1);
subunit_drive = exp(squeeze(sum(SU_covariates,1)));
for i_SU = 1:n_params
    del_lcif = zeros(1, length(SU_covariates));
    for i_loc = 1:121
        temp = pooling_weights(i_loc)*squeeze(pixels(i_SU,i_loc,:))'.*squeeze(subunit_drive(i_loc,:));
        del_lcif = del_lcif+temp;
    end
    del_lcif = conv(del_lcif, flip(time_filter), 'same');
    del_lcif = imresize(del_lcif, n_bins, 'nearest');
    del_LL(i_SU) = sum(del_lcif(spt))- dt * sum(exp(lcif).*del_lcif); 
end
del_penalty = zeros(n_params, 1);
for i_SU = 1:n_params
    del_penalty(i_SU) = f_penalty - penalty_strength*SU_params(i_SU)^2 + 2*penalty_strength*SU_params(i_SU);
end
grad = -del_LL + del_penalty;

% % Evaluate the hessian


%%

% Switch signs because using a minimizer  fmin
f       = -f_eval+f_penalty;
% grad    = -g_eval;
% Hess    = -H_eval;
% log_cif = lcif;
end
function  [f grad Hess log_cif] = glm_SU_optimizationfunction_LL(SU_params, SU_covariates, pooling_weights, spikebins, bin_duration, non_stim_lcif)
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
dt = bin_duration;
spt = spikebins;
n_bins = size(non_stim_lcif);
n_params = length(SU_params);

% Find Predicted Rate
pixels = SU_covariates;
for i = 1:n_params
    SU_covariates(i, :, :) = SU_covariates(i, :, :) * SU_params(i);
end
stim_lcif = log(pooling_weights'* exp(squeeze(sum(SU_covariates,1))));
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));

%% Evaluate the gradient
% gradient of the lcif
del_LL = zeros(n_params, 1);
subunit_drive = exp(squeeze(sum(SU_covariates,1)));
for i_SU = 1:n_params
    del_lcif = zeros(1, length(SU_covariates));
    for i_loc = 1:121
        top = pooling_weights(i_loc)*squeeze(pixels(i_SU,i_loc,:))'.*squeeze(subunit_drive(i_loc,:));
        del_lcif_top = del_lcif+top;
    end
    bottom = pooling_weights'*subunit_drive;
    del_lcif = del_lcif_top./bottom;
    del_lcif = imresize(del_lcif, n_bins, 'nearest');
    del_LL(i_SU) = sum(del_lcif(spt))- dt * sum(exp(lcif).*del_lcif); 
end

f = -f_eval;
grad = -del_LL;

% % Evaluate the hessian

end
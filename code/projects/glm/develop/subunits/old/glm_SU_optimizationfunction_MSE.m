function  [f grad Hess log_cif] = glm_SU_optimizationfunction(SU_params, SU_covariates, pooling_weights, spikebins, bin_duration, non_stim_lcif)
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
stim_lcif = pooling_weights'* exp(squeeze(sum(SU_covariates,1)));
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;
pred_rate = exp(lcif);

% Find Actual Rate
rec_rate = zeros(1,length(lcif));
rec_rate(spikebins) = 1;

% Evaluate the objective function (monotonic in log-likelihood)
f = sum((rec_rate - pred_rate).^2);
% 
% % Evaluate the gradient
% gradient of the lcif
del_MSE = zeros(n_params, 1);
subunit_drive = exp(squeeze(sum(SU_covariates,1)));
for i_SU = 1:n_params
    del_lcif = zeros(1, length(SU_covariates));
    for i_loc = 1:121
        temp = pooling_weights(i_loc)*squeeze(pixels(i_SU,i_loc,:))'.*squeeze(subunit_drive(i_loc,:));
        del_lcif = del_lcif+temp;
    end
    del_lcif = imresize(del_lcif, n_bins, 'nearest');
    del_MSE(i_SU) = sum(-2*exp(lcif).*(rec_rate - pred_rate).*del_lcif);
end
grad = del_MSE;

% % Evaluate the hessian

end
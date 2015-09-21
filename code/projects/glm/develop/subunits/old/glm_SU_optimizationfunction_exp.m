function  [f, grad, Hess, log_cif] = glm_SU_optimizationfunction_exp(SU_params, SU_covariates, pooling_weights, time_filter, spikebins, bin_duration, non_stim_lcif)
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
n_time = length(SU_covariates);
n_loc = length(pooling_weights);

%% Find Conditional Intensity and its log
pixels = SU_covariates;
subunit_drive = repmat(pooling_weights, [1, n_time]).*exp(squeeze(sum(repmat(SU_params, [1 n_loc n_time]).*SU_covariates)));
if ~(time_filter == 0)
    stim_lcif = conv(sum(subunit_drive), time_filter', 'full');
    stim_lcif = stim_lcif(1:n_time);
else
    stim_lcif = subunit_drive;
end
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;
clear sitm_lcif non_stim_lcif


% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));

%% Evaluate the gradient
% gradient of the lcif
del_lcif = squeeze(sum(pixels.*repmat(reshape(subunit_drive, [1 n_loc, n_time]), [n_params, 1,1]),2));
if ~(time_filter == 0)
    del_lcif = conv2(del_lcif, time_filter', 'full');
    del_lcif = del_lcif(1:n_time);
end
del_lcif = imresize(del_lcif, [n_params, n_bins(2)], 'nearest');
g_eval = sum(del_lcif(:,spt),2)'- dt * exp(lcif)*del_lcif';

%% Evaluate the hessian
H_eval = zeros(n_params);
for i = 1:n_params
    for j = 1:i
        di_dj_lcif = sum(squeeze(pixels(i,:,:).*pixels(j,:,:)).*subunit_drive);
        if ~(time_filter == 0)
            di_dj_lcif = conv(di_dj_lcif, time_filter, 'full');
            di_dj_lcif = di_dj_lcif(1:n_time);
        end
        di_dj_lcif = imresize(di_dj_lcif, n_bins, 'nearest');
        H_eval(i,j) = sum(di_dj_lcif(spt)) - dt * exp(lcif) * (di_dj_lcif + del_lcif(i,:).*del_lcif(j,:))';
        H_eval(j,i) = H_eval(i,j);
    end
end

%%
% ADD PENALTY
penalty_strength = 500;
%L2
% f_penalty = penalty_strength*(p'*p);
% del_penalty = 2*penalty_strength*p;
% saturating
f_penalty = penalty_strength * sum(atan(SU_params).^2);
del_penalty = 2*penalty_strength*atan(SU_params).*(SU_params.^2 + 1).^(-1);
%L1
% f_penalty = penalty_strength*sum(abs(p));
% del_penalty = penalty_strength*one(size(p));
% del_penalty(p<0) = -penalty_strength;

%%

% Switch signs because using a minimizer  fmin
f       = -f_eval+f_penalty;
grad    = -g_eval+del_penalty';
Hess    = -H_eval;
log_cif = lcif;
end

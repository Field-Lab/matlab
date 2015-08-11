function  [f, grad, Hess, log_cif] = glm_SU_time_optimizationfunction_exp(SU_time_params, SU_covariates, post_timefilter, pooling_weights, spikebins, bin_duration, non_stim_lcif)
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
tic 
%% Initialize
GLMPars = GLMParams; 
dt = bin_duration;
spt = spikebins;
n_bins = size(non_stim_lcif);
n_params = length(SU_time_params);
n_SU_params = GLMPars.subunit.size^2;
n_time_params = GLMPars.subunit.pretime_filter_frames;
n_time = length(SU_covariates);

SU_params = SU_time_params(1:n_SU_params);
pre_timefilter = SU_time_params((n_SU_params+1):end);
% penalty_strength = 0;


%% Find Conditional Intensity and its log
pixels = SU_covariates;
for i = 1:n_params
    SU_covariates(i, :, :) = SU_covariates(i, :, :) * SU_params(i);
end
subunit_drive = exp(conv(squeeze(sum(SU_covariates,1)), pre_timefilter, 'full'));
subunit_drive = subunit_drive(1:n_time);
if ~(post_timefilter == 0)
    stim_lcif = conv(pooling_weights'*subunit_drive, post_timefilter, 'full');
    stim_lcif = stim_lcif(1:n_time);
else 
    stim_lcif = pooling_weights'*subunit_drive;
end
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;
clear stim_lcif non_stim_lcif

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));

% Add penalty
% f_penalty = penalty_strength*SU_params'*SU_params;

%% Evaluate the gradient
% gradient of the lcif
SU_deriv = convn(pixels, pre_timefilter, 'full');
SU_deriv = SU_deriv(:, :, 1:n_time);

g_eval = zeros(n_params, 1);
del_lcif = zeros(n_params, n_bins(2));
for i_SU = 1:n_params
    temp_del_lcif = zeros(1,n_time);
    for i_loc = 1:121
        if i_SU > n_SU_params
            tau = i_SU - n_SU_params;
            T_deriv = [zeros(n_SU_params, 121, tau-1); SU_covariates(:,:,1:(end-tau)];
            temp = pooling_weights(i_loc)*(squeeze(SU_deriv(i_SU,i_loc,:)))'.*squeeze(subunit_drive(i_loc,:));
        else
            temp = pooling_weights(i_loc)*(squeeze(SU_deriv(i_SU,i_loc,:)))'.*squeeze(subunit_drive(i_loc,:));
        end
        temp_del_lcif = temp_del_lcif+temp;
    end
    if ~(time_filter == 0)
        temp_del_lcif = conv(temp_del_lcif, post_timefilter, 'full');
        temp_del_lcif = temp_del_lcif(1:n_time);
    end
    del_lcif(i_SU,:) = imresize(temp_del_lcif, n_bins, 'nearest');
    g_eval(i_SU) = sum(del_lcif(i_SU,spt))- dt * exp(lcif)*del_lcif(i_SU,:)';
end
% del_penalty = zeros(n_params, 1);
% for i_SU = 1:n_params
%     del_penalty(i_SU) = f_penalty - penalty_strength*SU_params(i_SU)^2 + 2*penalty_strength*SU_params(i_SU);
% end
%g_eval = del_LL;% - del_penalty;

%% Evaluate the hessian
H_eval = zeros(n_params);
for i = 1:n_params
    for j = 1:i
        di_dj_lcif = zeros(1,n_time);
        for i_loc = 1:121
            temp = pooling_weights(i_loc)*(squeeze(pixels(j,i_loc,:)).*squeeze(pixels(i,i_loc,:)))'.*squeeze(subunit_drive(i_loc,:));
            di_dj_lcif = di_dj_lcif+temp;
        end
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

% Switch signs because using a minimizer  fmin
f       = -f_eval;%+f_penalty;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end
function  [f, grad, Hess, log_cif] = glm_SU_time_optimizationfunction_exp(SU_time_params, SU_covariates, pooling_weights,post_timefilter, spikebins, bin_duration, non_stim_lcif)
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
n_loc = length(pooling_weights);

SU_params = SU_time_params(1:n_SU_params);
pre_timefilter = SU_time_params((n_SU_params+1):end);
post_timefilter = 0;
% penalty_strength = 0;


%% Find Conditional Intensity and its log

pixels = SU_covariates;

% timefilt_pixels = T*S
timefilt_pixels = convn(SU_covariates, reshape(pre_timefilter, [1, 1, n_time_params]), 'full');
timefilt_pixels = timefilt_pixels(:,:,1:n_time);

% SU_cov = K*S
SU_covariates = squeeze(sum(SU_covariates.*repmat(SU_params, [1 n_loc n_time]),1));

% SU_drive = exp(K*S)
subunit_drive = repmat(pooling_weights, [1, n_time]).*exp(squeeze(sum(repmat(SU_params, [1 n_loc n_time]).*timefilt_pixels)));

% Convolve with post filter if necessary
if ~(post_timefilter == 0)
    stim_lcif = conv(sum(subunit_drive), post_timefilter, 'full'); % will need to check dim on this
    stim_lcif = stim_lcif(1:n_time);
else 
    stim_lcif = sum(subunit_drive);
end
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;
clear stim_lcif non_stim_lcif

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));

% Add penalty
% f_penalty = penalty_strength*SU_params'*SU_params;

%% Evaluate the gradient
% gradient of the lcif
del_lcif = zeros(n_params, n_time);
del_lcif(1:n_SU_params,:) = squeeze(sum(timefilt_pixels.*repmat(reshape(subunit_drive, [1 n_loc, n_time]), [n_SU_params, 1,1]),2));
for tau = 1:n_time_params
    del_lcif(n_SU_params+tau,:) = sum([zeros(n_loc, tau-1) SU_covariates(:,tau:end)].*subunit_drive);
end
if ~(post_timefilter == 0)
    del_lcif = convn(del_lcif, reshape(post_timefilter, [1, 30]), 'full');
    del_lcif = del_lcif(:,1:n_time);
end
del_lcif = imresize(del_lcif, [n_params n_bins(2)], 'nearest');
g_eval = sum(del_lcif(:,spt),2)'- dt * exp(lcif)*del_lcif';

% del_penalty = zeros(n_params, 1);
% for i_SU = 1:n_params
%     del_penalty(i_SU) = f_penalty - penalty_strength*SU_params(i_SU)^2 + 2*penalty_strength*SU_params(i_SU);
% end
%g_eval = del_LL;% - del_penalty;

%% Evaluate the hessian
% NEED TO FIX THIS
H_eval = zeros(n_params);
for i = 1:n_params
    for j = 1:i
        
        % Both spatial
        if i <= n_SU_params
            di_dj_lcif = sum(squeeze(timefilt_pixels(i,:,:).*timefilt_pixels(j,:,:)).*subunit_drive);
        elseif i > n_SU_params
            tau1 = i-n_SU_params;
            
            % both temporal
            if j > n_SU_params
                tau2 = j-n_SU_params;
                T_deriv = [zeros(n_loc, tau1-1) SU_covariates(:,tau1:end)].*[zeros(n_loc, tau2-1) SU_covariates(:,tau2:end)];
                di_dj_lcif = sum(T_deriv.*subunit_drive);
                
            else% mixed
                T_deriv = [zeros(n_loc, tau1-1) SU_covariates(:,tau1:end)];
                di_dj_lcif = sum(squeeze(timefilt_pixels(j,:,:)).*T_deriv.*subunit_drive) +sum(subunit_drive.*[zeros(n_loc, tau1-1) squeeze(pixels(j,:,tau1:end))]);
            end
        end
        
        if ~(post_timefilter == 0)
            di_dj_lcif = conv(di_dj_lcif, post_timefilter, 'full');
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
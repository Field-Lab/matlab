function  [f grad Hess log_cif] = glm_SU_time_optimizationfunction(SU_time_params, SU_covariates, pooling_weights, spikebins, bin_duration, non_stim_lcif)
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
dt = bin_duration;
spt = spikebins;
n_bins = size(non_stim_lcif);
n_params = length(SU_time_params);
n_SU_params = 8;
n_time_params = 30;
n_time = length(SU_covariates);

SU_params = SU_time_params(1:n_SU_params);
time_filter = SU_time_params((n_SU_params+1):end);
% penalty_strength = 0;

%% Find Conditional Intensity and its log
pixels = SU_covariates(1:n_SU_params,:,:)-repmat(SU_covariates(n_SU_params+1,:,:),[n_SU_params 1 1]); % for sum 0 constraint

for i = 1:n_SU_params
    SU_covariates(i, :, :) = SU_covariates(i, :, :) * SU_params(i);
end
SU_covariates(n_SU_params+1,:,:) = SU_covariates(n_SU_params+1,:,:) * -sum(SU_params); % the last SU value is set so the sum is 0 
subunit_drive = exp(conv2(squeeze(sum(SU_covariates,1)),flip(time_filter)', 'full')); % sum_tau T_tau * K * S_loc
subunit_drive = subunit_drive(:,1:n_time);

stim_lcif = pooling_weights'* subunit_drive;
% stim_lcif = conv(stim_lcif, flip(time_filter), 'same'); % put this back
% in for second time filter
lcif = imresize(stim_lcif, n_bins, 'nearest') + non_stim_lcif;

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));
toc
% Add penalty
% f_penalty = penalty_strength*SU_params'*SU_params;

%% Evaluate the gradient
% gradient of the lcif
g_eval = zeros(n_params, 1);
del_lcif = zeros(n_params, n_bins(2));
for i_SU = 1:n_SU_params
    temp_del_lcif = zeros(1,length(pixels));
    for i_loc = 1:length(pooling_weights)
        temp_conv = conv(squeeze(pixels(i_SU,i_loc,:)), time_filter, 'full');
        temp = pooling_weights(i_loc)*temp_conv(1:n_time)'.*squeeze(subunit_drive(i_loc,:));
        temp_del_lcif = temp_del_lcif+temp;
    end
    % temp_del_lcif = conv(temp_del_lcif, flip(time_filter), 'same');
    del_lcif(i_SU,:) = imresize(temp_del_lcif, n_bins, 'nearest');
    g_eval(i_SU) = sum(del_lcif(i_SU,spt))- dt * exp(lcif)*del_lcif(i_SU,:)';
end

SU_cov_shift = zeros([n_time_params,size(SU_covariates)]);
for i_time = 1:n_time_params
    % THIS IS TOO SLOW: NEED TO CHANGE
    time_shift = cat(2, zeros(length(pooling_weights), i_time), squeeze(sum(SU_covariates(:,:,1:(end-i_time)))));
    temp_del_lcif = zeros(1,length(pixels));
    for i_loc = 1:length(pooling_weights)
        temp = pooling_weights(i_loc)*time_shift(i_loc,:).*squeeze(subunit_drive(i_loc,:));
        temp_del_lcif = temp_del_lcif+temp;
    end
    % temp_del_lcif = conv(temp_del_lcif, flip(time_filter), 'full');
    % temp_del_lcif = temp_del_lcif(1:n_time);
    del_lcif(i_SU,:) = imresize(temp_del_lcif, n_bins, 'nearest');
    g_eval(i_SU) = sum(del_lcif(i_SU,spt))- dt * exp(lcif)*del_lcif(i_SU,:)';
end
% del_penalty = zeros(n_params, 1);
% for i_SU = 1:n_params
%     del_penalty(i_SU) = f_penalty - penalty_strength*SU_params(i_SU)^2 + 2*penalty_strength*SU_params(i_SU);
% end
%g_eval = del_LL;% - del_penalty;

toc

%% Evaluate the hessian
H_eval = zeros(n_params);
for i = 1:n_params
    if i > n_SU_params
        tau = i-n_SU_params;
        time_shift_i = cat(2, zeros(length(pooling_weights), tau), squeeze(sum(SU_covariates(:,:,1:(end-tau)))));
    end
    for j = 1:i
        if j > n_SU_params
            tau = j-n_SU_params;
            time_shift_j = cat(2, zeros(length(pooling_weights), tau), squeeze(sum(SU_covariates(:,:,1:(end-tau)))));
        end
        di_dj_lcif = zeros(1,n_time);
        % both spatial params
        if i <= n_SU_params && j <= n_SU_params
            for i_loc = 1:length(pooling_weights)
                temp_conv = conv(squeeze(pixels(i,i_loc,:)), time_filter, 'full').*conv(squeeze(pixels(j,i_loc,:)), time_filter, 'full');
                temp = pooling_weights(i_loc)*temp_conv(1:n_time)'.*squeeze(subunit_drive(i_loc,:));
                di_dj_lcif = di_dj_lcif+temp;
            end
            % one spatial j and one temporal i params
        elseif i > n_SU_params && j <= n_SU_params
            for i_loc = 1:length(pooling_weights)
                temp_conv = conv(squeeze(pixels(j,i_loc,:)), time_filter, 'full');
                temp1 = pooling_weights(i_loc)*(time_shift_i(i_loc,:).*temp_conv(1:n_time)').*squeeze(subunit_drive(i_loc,:));
                temp2 = pooling_weights(i_loc)*squeeze(pixels(j,i_loc,:))'.*squeeze(subunit_drive(i_loc,:));
                di_dj_lcif = di_dj_lcif+temp1+temp2;
            end
            % both temporal params
        else
            for i_loc = 1:length(pooling_weights)
                temp_conv = time_shift_i(i_loc,:).*time_shift_j(i_loc,:);
                temp = pooling_weights(i_loc)*temp_conv.*squeeze(subunit_drive(i_loc,:));
                di_dj_lcif = di_dj_lcif+temp;
            end
        end
        
        % di_dj_lcif = conv(di_dj_lcif, flip(time_filter), 'same');
        di_dj_lcif = imresize(di_dj_lcif, n_bins, 'nearest');
        H_eval(i,j) = sum(di_dj_lcif(spt)) - dt * exp(lcif) * (di_dj_lcif + del_lcif(i,:).*del_lcif(j,:))';
        H_eval(j,i) = H_eval(i,j);
    end
end
toc

%%

% Switch signs because using a minimizer  fmin
f       = -f_eval;%+f_penalty;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end
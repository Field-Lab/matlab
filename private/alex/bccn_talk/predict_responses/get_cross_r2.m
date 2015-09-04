function [r2, norm_asr1, norm_asr2] = get_cross_r2(datarun, cellID, trig_threshold, n_repeats, sta_params, norm_method, refresh)

mean_asr = get_mean_spike_rate(datarun{1}, cellID, trig_threshold, n_repeats, refresh(1));
mean_asr = mean_asr(sta_params.length+1:end);
eval(['norm_asr1 = mean_asr/', norm_method,'(mean_asr);']);

% get actual response: datarun 2
mean_asr = get_mean_spike_rate(datarun{2}, cellID, trig_threshold, n_repeats, refresh(2));
mean_asr = mean_asr(sta_params.length+1:end);
eval(['norm_asr2 = mean_asr/', norm_method,'(mean_asr);']);

% r2: datarun 1 predicted by datarun 2
sse = sum((norm_asr2-norm_asr1).^2); % variance between model and data
sst = sum((norm_asr1-mean(norm_asr1)).^2); % variance in the data
r2 = 1 - sse/sst; % fraction of data variance exlpained by model

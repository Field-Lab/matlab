function [filtered_lum,filter_log,times,luminence_movie] = luminence_signal_bank(maskedMovdd_sliced,fitGMLM,WN_data)

tau_list = [20]/1000;% [20:50:520]/1000; % in, 1 ms 300/1000;

times = [0:max(tau_list)*120*3.5]*1/120; % in seconds

filter_log = zeros(length(times),length(tau_list));

for itau = 1:length(tau_list)
    tau = tau_list(itau);
filter_tau = exp(-times/tau);
filter_log(:,itau) = filter_tau;

end

%filter_log = filter_log-repmat(mean(filter_log,1),[size(filter_log,1),1]);
figure;
plot(times,filter_log);
title('filters');

[q,r] = qr(filter_log,0);
filter_log = q;
figure;
plot(times,filter_log);
title('filter basis');

filtered_lum = []%;zeros(length(tau_list),size(maskedMovdd_sliced,2));

%% mean in RF
%  luminence_movie = [mean(maskedMovdd_sliced,1)];
%  filtered_lum = filter_lum(luminence_movie,filter_log);
 
%% mean change in RF 
%   luminence_movie =[0, diff(mean(maskedMovdd_sliced,1))];
% filtered_lum2 = filter_lum(luminence_movie,filter_log);
% filtered_lum = [filtered_lum;filtered_lum2];

%% mean in SUs 

nSU = length(fitGMLM.Linear.filter);

for isu=1:nSU
luminence_movie = fitGMLM.Linear.filter{isu}'*maskedMovdd_sliced;
filtered_lum3 = filter_lum(luminence_movie,filter_log);
filtered_lum = [filtered_lum ; filtered_lum3];
end

%% mean change in SUs


nSU = length(fitGMLM.Linear.filter);

for isu=1:nSU
luminence_movie = fitGMLM.Linear.filter{isu}'*maskedMovdd_sliced;
luminence_movie = [ 0 , diff(luminence_movie,1)];
filtered_lum3 = filter_lum(luminence_movie,filter_log);
filtered_lum = [filtered_lum ; filtered_lum3];
end

end
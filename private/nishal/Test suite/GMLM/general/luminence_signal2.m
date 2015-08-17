function [filtered_lum,filter_log,times,luminence_movie] = luminence_signal2(maskedMovdd_sliced)

tau_list =[20]/1000; % in, 1 ms 

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

filtered_lum = zeros(length(tau_list),size(maskedMovdd_sliced,2));
%luminence_movie =[0,diff(mean(maskedMovdd_sliced,1))]./((mean((maskedMovdd_sliced),1))-min((maskedMovdd_sliced(:))) +0.01);
% luminence_movie =[0,diff(mean(maskedMovdd_sliced,1))];
luminence_movie = [mean(maskedMovdd_sliced,1)];

for itau = 1: size(filter_log,2)
    filter_tau = filter_log(:,itau);
luminence_movie_ex = zeros(length(luminence_movie) + length(filter_tau)-1,1);
luminence_movie_ex(length(filter_tau):end) = luminence_movie;

filtered_lum(itau,:) = conv(luminence_movie_ex,filter_tau,'valid');
end

% 
% filtered_lum2 = zeros(length(tau_list),size(maskedMovdd_sliced,2));
% luminence_movie =[0, diff(mean(maskedMovdd_sliced,1))];
% for itau = 1: size(filter_log,2)
%     filter_tau = filter_log(:,itau);
% luminence_movie_ex = zeros(length(luminence_movie) + length(filter_tau)-1,1);
% luminence_movie_ex(length(filter_tau):end) = luminence_movie;
% 
% filtered_lum2(itau,:) = conv(luminence_movie_ex,filter_tau,'valid');
% end
% %filtered_lum = filtered_lum2;
% 
% %filtered_lum = [filtered_lum;filtered_lum2];

end
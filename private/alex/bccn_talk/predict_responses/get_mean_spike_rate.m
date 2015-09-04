function [mean_sr, inter_corr] = get_mean_spike_rate(spikes, triggers, trig_threshold, n_repeats, fps, meth)


trial_begins=round([0; triggers(find(diff(triggers)>trig_threshold)+1)]*fps);

spikes=ceil(spikes*fps);
spikes(spikes>max(diff(trial_begins))*n_repeats)=[];
spikes(spikes<1) = [];

% get continuous spike rate 
spike_rate=zeros(max(diff(trial_begins))*n_repeats,1);
while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
    spikes(ia)=[];
end
    
if  strcmp(meth, 'smooth') % convolve with a kernel
    sr=120;
    sig=40;
    st=10000/sr*6.2/(60*sig);
    time_list=-3.1:st:3.1;
    kern=zeros(1,length(time_list));
    for i=1:length(time_list)
        kern(i)=250/sig*exp((1-time_list(i)^2)/2);
    end
    spike_rate = conv(spike_rate,kern,'same');
end

% get mean spike rate across repeats
mean_sr=0;
for i=1:n_repeats
    mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+min(diff(trial_begins)));
end
mean_sr=mean_sr/n_repeats;

% get half repeats correaltion
tmp = 0;
tmp1 = 0;
for i=1:2:n_repeats-1
    tmp = tmp + spike_rate(trial_begins(i)+1:trial_begins(i)+min(diff(trial_begins)));
    tmp1 = tmp1 + spike_rate(trial_begins(i+1)+1:trial_begins(i+1)+min(diff(trial_begins)));
end
inter_corr = corr(tmp, tmp1);

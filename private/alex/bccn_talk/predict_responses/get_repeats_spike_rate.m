function repeats = get_repeats_spike_rate(spikes, triggers, trig_threshold, n_repeats, fps)


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
    
% get mean spike rate across repeats
repeats=zeros(n_repeats, min(diff(trial_begins)));
for i=1:n_repeats
    repeats(i,:) = spike_rate(trial_begins(i)+1:trial_begins(i)+min(diff(trial_begins)));
end



function spikes_adj = align_spikes_triggers(spikes, triggers, frames_per_trigger, monitor_refresh)
    %spikes_adj=spikes;
    spikes_adj = [];
    n_block=0;
    block_length = frames_per_trigger/monitor_refresh;
    for i=1:(length(triggers))
        actual_t_start=triggers(i);
        supposed_t_start=n_block*block_length;
        idx1=spikes > actual_t_start;
        idx2=spikes < actual_t_start+block_length;%triggers(i+1);
        %spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
        spikes_adj=[spikes_adj; spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start];
        n_block=n_block+1;
    end
end
function spikes_adj = align_spikes_triggers(spikes, triggers, frames_per_trigger, monitor_refresh)
    spikes_adj=spikes;
    n_block=0;
    for i=1:(length(triggers)-1)
        actual_t_start=triggers(i);
        supposed_t_start=n_block*frames_per_trigger/monitor_refresh;
        idx1=spikes > actual_t_start;
        idx2=spikes < triggers(i+1);
        spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
        n_block=n_block+1;
    end
end

function unblocked_raster(datarun, cell, start_times, raster_length)
    spikes=datarun.spikes{cell};
    n_raster=length(start_times);
    hFig1=figure;
    set(hFig1, 'Position', [100 100 800 250])
    hold on
    for i=1:n_raster
        start=start_times(i);
        spikes_temp=spikes(spikes > start & spikes < (start+raster_length));
        plot(spikes_temp-start,i*ones(length(spikes_temp),1),'b.');
    end
    ylim([1 n_raster])
    ylabel('Trials')
    xlabel('Seconds')
    
end
function slope = amacrine_wave_plot(prepped_data, datarun_class)

cids = get_cell_indices(datarun_class, prepped_data.p.Results.cell_spec);

figure; hold on
slope = zeros(20, 2);
for i_saccade = 0:19
    for i_cell = 1:length(cids)
        spikes = cell2mat(prepped_data.testspikes(:,i_cell));
        centers = flip(datarun_class.vision.sta_fits{cids(i_cell)}.mean);
        hold on; subplot(2,1,1); plot(spikes, centers(1)*ones(length(spikes), 1), '.k')
        title(['Saccade ' num2str(i_saccade)])
        xlim([i_saccade i_saccade+0.25])
        hold on; subplot(2,1,2); plot(spikes, centers(2)*ones(length(spikes), 1), '.k')
        xlim([i_saccade i_saccade+0.25])
    end
    pause()
end

end
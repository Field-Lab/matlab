
cids = get_cell_indices(datarun_class, 'Off Amacrine');

hFig = figure; hold on
set(hFig, 'PaperOrientation', 'landscape')
set(hFig, 'Position', [100 100 900 250])
% slope = zeros(20, 2);

for i_cell = 1:length(cids)
    spikes = cell2mat(prepped_data.testspikes(:,i_cell));
    spike_frame = floor(spikes*120);
    for i = 30*120
        spike_rate(i_cell,i)=sum(spike_frame == i);
    end
    centers = flip(datarun_class.vision.sta_fits{cids(i_cell)}.mean);
    %hold on; subplot(2,1,1); plot(spike_rate(i_cell,:) + 10/scale*centers(1), 'k')
    %hold on; subplot(2,1,2); plot(spike_rate(i_cell,:) + 10/scale*centers(2), 'k')
%     plotting_matrix1 = repmat(spikes, [2,1]); 
%     plotting_matrix2 = [10/scale*centers(1)*ones(length(spikes),1)+10; 10/scale*centers(1)*ones(length(spikes),1)-10];
%     plotting_matrix = [plotting_matrix1 plotting_matrix2];
%     hold on; subplot(2,1,1); plot(plotting_matrix', 'k');
    hold on; subplot(2,1,1); plot(spikes,10*centers(1)*ones(length(spikes),1), 'k.')
    hold on; subplot(2,1,2); plot(spikes,10*centers(2)*ones(length(spikes),1), 'k.')
end

for i_saccade = 0:29
    %snippet =  spike_rate(:, (1+i_saccade*bin):(i_saccade*bin+40*scale));
    %[~,i] = max(diff(snippet(:,1:8)'));
%     for i_cell = 1:length(cids)
%         hold on; subplot(2,1,1); plot(snippet(i_cell, :) + 10*centers(i_cell,1),'k');
%         %hold on; subplot(2,1,1); plot(10/scale*centers(i_cell,1),i(i_cell),  'r.')
%         hold on; subplot(2,1,2); plot(snippet(i_cell, :) + 10*centers(i_cell,2),'k');
%         %hold on; subplot(2,1,2); plot(10/scale*centers(i_cell,2),i(i_cell), 'r.')
%     end
    subplot(2,1,1); xlim([i_saccade i_saccade+0.25])
    title(['Saccade ' num2str(i_saccade) 'Exp 05'])
    subplot(2,1,2); xlim([i_saccade i_saccade+0.25])
    print(['/Users/Nora/Desktop/saccade_images/saccade10_' num2str(i_saccade)], '-dpdf')
end

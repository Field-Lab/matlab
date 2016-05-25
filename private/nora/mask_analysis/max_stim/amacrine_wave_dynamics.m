on = [0 0];
off = [0 0];
large = [0 0];


%%

%{
datarun = load_data('2015-10-29-2/map-from-data036-stream/data040/data040'); % masked 
datarun_LES = load_data('2015-10-29-2/map-from-data036-stream/data041/data041'); % LES 
cells = [1801,3287,4907,6393,6841,7263];
type = {'off','off', 'on', 'on','off','on'};
%}
%{
datarun = load_data('2015-10-29-2/map-from-data036-stream/data042/data042'); % masked 
datarun_LES = load_data('2015-10-29-2/map-from-data036-stream/data043/data043'); % LES 
cells = [753,2134,2687,3677,4246,5222,6707];
type = {'on','on','on','off','off','on','off'};
%}
%{
datarun = load_data('2015-10-29-7/mapped-from-data001-streamed/data004/data004'); % masked 
datarun_LES = load_data('2015-10-29-7/mapped-from-data001-streamed/data005/data005'); % LES 
cells = [1411,2072,3152,4397,6181,7217];
type = {'on','off','on','off','on','off'};
%}
%%{
datarun = load_data('2015-10-29-7/mapped-from-data001-streamed/data006/data006'); % masked 
datarun_LES = load_data('2015-10-29-7/mapped-from-data001-streamed/data007/data007'); % LES 
cells = [1727,2120,3616,4611,5778,7608];
type = {'off','large','off','large','off','on'};
%}

datarun = load_neurons(datarun);
datarun_LES = load_neurons(datarun_LES);
mask = interleaved_data_prep(datarun, 1800, 40, 'cell_spec', cells, 'visual_check', 0);
LES = interleaved_data_prep(datarun_LES, 1800, 40, 'cell_spec', cells, 'visual_check', 0);

bins_per_frame = 1;
bins = 1800*bins_per_frame;

PSTH_LES=zeros(bins,6);
PSTH_mask=zeros(bins,6);
for i_cell = 1:6
    spike_frames_LES = floor(cell2mat(LES.testspikes(:,i_cell))*120*bins_per_frame);
    spike_frames_mask = floor(cell2mat(mask.testspikes(:,i_cell))*120*bins_per_frame);
    for i = 1:(1800*bins_per_frame)
        PSTH_LES(i,i_cell) = sum(spike_frames_LES == i);
        PSTH_mask(i,i_cell) = sum(spike_frames_mask == i);
    end
end
PSTH_LES = conv2(PSTH_LES, gausswin(10), 'same');
PSTH_mask = conv2(PSTH_mask, gausswin(10), 'same');

%%
for i = 1:100
    load(['/Volumes/Lab/Users/Nora/new_stim_nora/NSEM intervals/matfiles/movie_chunk_' num2str(i) '.mat'])
    interval(i) = size(movie_chunk, 3);
end
default_colors = get(gca,'ColorOrder');
interval_time = cumsum(interval/120);



%%
hFig1=figure;
default_colors = get(gca,'ColorOrder');
set(hFig1, 'Position', [100 100 1800 250])
for i_cell =1:6
    plot(PSTH_mask(:,i_cell))%,PSTH_LES(:,i), '.');
    hold on
    plot(PSTH_LES(:,i_cell));
    % plot([0 150], [0 150]);
    i = 1;
    time = interval_time(1);
    while time < 15
        plot([120*time 120*time], [0 200], 'Color', default_colors(3,:))
        i = i+1;
        time = interval_time(i);
    end
    hold off
    ylim([0 200])
    title([type{i_cell} num2str(cells(i_cell))])
    pause()
end

%% spikes per image
frame_changes = cumsum(interval);
for i_cell = 1:6
    i = 1;
    time = 0;
    while time < 15
        try
            mask_spike_count(i) = max(PSTH_mask(frame_changes(i):frame_changes(i+1), i_cell));
            % mask_spike_count_error(i) = std(PSTH_mask(frame_changes(i):frame_changes(i+1), i_cell));
            LES_spike_count(i) = max(PSTH_LES(frame_changes(i):frame_changes(i+1), i_cell));
            % LES_spike_count_error(i) = std(PSTH_LES(frame_changes(i):frame_changes(i+1), i_cell));
        catch
            mask_spike_count(i) = max(PSTH_mask(frame_changes(i):end, i_cell));
            % mask_spike_count_error(i) = std(PSTH_mask(frame_changes(i):end, i_cell));
            LES_spike_count(i) = max(PSTH_LES(frame_changes(i):end, i_cell));
            % LES_spike_count_error(i) = std(PSTH_LES(frame_changes(i):end, i_cell));
        end
        i = i+1;
        time = frame_changes(i)/120;
    end 
    switch type{i_cell}
        case 'on'
            on = [on; mask_spike_count', LES_spike_count'];%, mask_spike_count_error', LES_spike_count_error'];
        case 'off'
            off = [off; mask_spike_count', LES_spike_count'];%, mask_spike_count_error', LES_spike_count_error'];
        case 'large'
            large = [large; mask_spike_count', LES_spike_count'];%, mask_spike_count_error', LES_spike_count_error'];  
    end
    %errorbarxy(mask_spike_count, LES_spike_count, mask_spike_count_error/sqrt(40),LES_spike_count_error/sqrt(40))
%     plot(mask_spike_count, LES_spike_count, '.');
%     hold on
%     plot([0 1200], [0 1200])
%     hold off
%     axis equal
%     xlim([0 1200])
%     ylim([0 1200])
%     title([type{i_cell} num2str(cells(i_cell))])
%     xlabel('Original Stimulus');
%     ylabel('Linear Equivalent Stimulus');
%     pause()
end


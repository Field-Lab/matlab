on = [0 0];
off = [0 0];
large = [0 0];

%{
for i = 1:100
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i) '.mat'])
    interval(i) = size(movie_chunk, 3);
end
default_colors = get(gca,'ColorOrder');
interval_time = cumsum(interval/120);
%}

%%
% 1 3-5 
for exp = [1 3:5]
    
    switch exp
        case 1
            datarun = load_data('2015-10-29-2/map-from-data036-stream/data040/data040'); % masked
            datarun_LES = load_data('2015-10-29-2/map-from-data036-stream/data041/data041'); % LES
            datarun_full = load_data('2015-10-29-2/map-from-data036-stream/data039/data039'); % full
            cells = [1801,3287,4907,6393,6841,7263];
            type = {'off','off', 'on', 'on','off','on'};
            attributes = {'29-1', };
        case 2
            datarun = load_data('2015-10-29-2/map-from-data036-stream/data042/data042'); % masked
            datarun_LES = load_data('2015-10-29-2/map-from-data036-stream/data043/data043'); % LES
            datarun_full = load_data('2015-10-29-2/map-from-data036-stream/data039/data039'); % full
            cells = [753,2134,2687,3677,4246,5222,6707];
            type = {'on','on','on','off','off','on','off'};
        case 3
            datarun = load_data('2015-10-29-7/mapped-from-data001-streamed/data004/data004'); % masked
            datarun_LES = load_data('2015-10-29-7/mapped-from-data001-streamed/data005/data005'); % LES
            datarun_class = load_data('2015-10-29-7/mapped-from-data001-streamed/data003/data003'); %full
            cells = [1411,2072,3152,4397,6181,7217];
            type = {'on','off','on','off','on','off'};
        case 4
            datarun = load_data('2015-10-29-7/mapped-from-data001-streamed/data006/data006'); % masked
            datarun_LES = load_data('2015-10-29-7/mapped-from-data001-streamed/data007/data007'); % LES
            datarun_class = load_data('2015-10-29-7/mapped-from-data001-streamed/data003/data003'); %full
            cells = [1727,2120,3616,4611,5778,7608];
            type = {'off','large','off','large','off','on'};
        case 5
            datarun = load_data('2015-11-09-1/map-from-data000-stream/data003/data003'); % masked
            datarun_LES = load_data('2015-11-09-1/map-from-data000-stream/data004/data004'); % LES
            datarun_full = load_data('2015-11-09-1/map-from-data000-stream/data002/data002'); %full
            cells = [182,751,1684, 1756, 2313, 3333, 4201, 4367, 6272, 7276];
            type = {'on','off','on','off','on', 'off','on','off', 'on', 'on'};
        case 6
            datarun = load_data('2015-11-09-1/map-from-data000-stream/data006/data006'); % masked
            datarun_LES = load_data('2015-11-09-1/map-from-data000-stream/data007/data007'); % LES
            datarun_full = load_data('2015-11-09-1/map-from-data000-stream/data008/data008'); %full
            cells = [182,751,1684, 1756, 2313, 3333, 4201, 4367, 6272, 7276];
            type = {'on','off','on','off','on', 'off','on','off', 'on', 'on'};
        case 7
            datarun = load_data('2015-11-09-8/map-from-data006-stream/data008/data008'); % masked
            datarun_LES = load_data('2015-11-09-8/map-from-data006-stream/data009/data009'); % LES
            datarun_full = load_data('2015-11-09-8/map-from-data006-stream/data010/data010'); % LES
            cells = [3,1147,1188,2208,2643,2866,3512,4863,5060,5915,6706,7636];
            type = {'off','on','on','on','on', 'off','off','off','on','on', 'off', 'off'};
        case 8
            datarun = load_data('2015-11-09-8/map-from-data006-stream/data012/data012'); % masked
            datarun_LES = load_data('2015-11-09-8/map-from-data006-stream/data013/data013'); % LES
            datarun_full = load_data('2015-11-09-8/map-from-data006-stream/data014/data014'); % LES
            cells = [3,1147,1188,2208,2643,2866,3512,4863,5060,5915,6706,7636];
            type = {'off','on','on','on','on', 'off','off','off','on','on', 'off', 'off'};
    end
    datarun = load_neurons(datarun);
    datarun_LES = load_neurons(datarun_LES);
    mask = interleaved_data_prep(datarun, 1800, 40, 'cell_spec', cells, 'visual_check', 0);
    LES = interleaved_data_prep(datarun_LES, 1800, 40, 'cell_spec', cells, 'visual_check', 0);
    
    bins_per_frame = 1;
    bins = 1800*bins_per_frame;
    
    n_cells = length(type);
    PSTH_LES=zeros(bins,n_cells);
    PSTH_mask=zeros(bins,n_cells);
    for i_cell = 1:n_cells
        spike_frames_LES = floor(cell2mat(LES.testspikes(:,i_cell))*120*bins_per_frame);
        spike_frames_mask = floor(cell2mat(mask.testspikes(:,i_cell))*120*bins_per_frame);
        for i = 1:(1800*bins_per_frame)
            PSTH_LES(i,i_cell) = sum(spike_frames_LES == i);
            PSTH_mask(i,i_cell) = sum(spike_frames_mask == i);
        end
        %     plot(PSTH_LES(:,i_cell))
        %     hold on; plot(PSTH_mask(:,i_cell))
        %     hold off
        %     pause()
        if (exp ==4 && (i_cell == 3 || i_cell == 2)) || (exp ==5 && i_cell == 3)
            IDP_plot_raster(mask, i_cell); title(['mask ' num2str(i_cell) num2str(exp) type{i_cell}])
            ylim([ 1 30])
            set(gcf, 'Position', [1 1 12 2])
            axis off
            exportfig(gcf, ['/Users/Nora/Desktop/' type{i_cell} '_spot_raster.eps'], 'Bounds', 'loose', 'Color', 'rgb');
            IDP_plot_raster(LES, i_cell);title('LES')
            ylim([ 1 30])
            set(gcf, 'Position', [1 1 12 2])
            axis off
            exportfig(gcf, ['/Users/Nora/Desktop/' type{i_cell} '_LES_raster.eps'], 'Bounds', 'loose', 'Color', 'rgb');
            close all
        end
    end
    PSTH_LES = conv2(PSTH_LES, gausswin(10), 'same');
    PSTH_mask = conv2(PSTH_mask, gausswin(10), 'same');
    %{
%%
hFig1=figure;
set(hFig1, 'Position', [100 100 1800 250])
for i =1:n_cells
    plot(PSTH_LES(:,i))
    hold on
    plot(PSTH_mask(:,i))
    hold off
    title(type{i})
    legend('Full screen', 'Center Only')
    pause()
end



for i_cell = 1:n_cells
    hold on
    for i = 1:40
        plot(mask.testspikes{i}, i*ones(length(mask.testspikes{i})), '.k')
    end
    hold off
    pause()
end
    %}
    
    % %%
    % hFig1=figure;
    % default_colors = get(gca,'ColorOrder');
    % set(hFig1, 'Position', [100 100 1800 250])
    % for i_cell =1:n_cells
    %     plot(PSTH_mask(:,i_cell))%,PSTH_LES(:,i), '.');
    %     hold on
    %     plot(PSTH_LES(:,i_cell));
    %     % plot([0 150], [0 150]);
    %     i = 1;
    %     time = interval_time(1);
    %     while time < 15
    %         plot([120*time 120*time], [0 200], 'Color', default_colors(3,:))
    %         i = i+1;
    %         time = interval_time(i);
    %     end
    %     hold off
    %     ylim([0 200])
    %     title([type{i_cell} num2str(cells(i_cell))])
    %     pause()
    % end
    
    %
    % spikes per image
    frame_changes = cumsum(interval);
    for i_cell = 1:n_cells
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
end

%%
on = log(on);
off = log(off);
large = log(large);
frame_size = 5;
figure;
plot(on(:,1), on(:,2), '.', 'MarkerSize', 10)
set(gcf, 'Position', [1 1 4 4])
axis equal
hold on; plot([0, frame_size], [0, frame_size])
xlim([0 frame_size])
ylim([0 frame_size])
title('ON')
exportfig(gcf, ['/Users/Nora/Desktop/' type{i_cell} '_ON_counts.eps'], 'Bounds', 'loose', 'Color', 'rgb');


figure;
plot(off(:,1), off(:,2), '.', 'MarkerSize', 10)
set(gcf, 'Position', [1 1 4 4])
axis equal
hold on; plot([0, frame_size], [0, frame_size])
xlim([0 frame_size])
ylim([0 frame_size])
title('OFF')
exportfig(gcf, ['/Users/Nora/Desktop/' type{i_cell} '_OFF_counts.eps'], 'Bounds', 'loose', 'Color', 'rgb');

figure;
plot(large(:,1), large(:,2), '.', 'MarkerSize', 10)
set(gcf, 'Position', [1 1 4 4])
axis equal
hold on; plot([0, frame_size], [0, frame_size])
xlim([0 frame_size])
ylim([0 frame_size])
title('ON Large')
exportfig(gcf, ['/Users/Nora/Desktop/' type{i_cell} '_LARGE_counts.eps'], 'Bounds', 'loose', 'Color', 'rgb');




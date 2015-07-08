% runGenSignal

clear all
close all

run_opts.date='2010-09-24-0/';
run_opts.concatname='data001-nwpca-cr';
% run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];

run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-80x60.xml';

run_opts.num_frames = 19;
run_opts.num_bins = 10; % Number of bins to use for the generator vs spikes graph

% cell_specification_large = [21]; % visionID
cell_specification = {'OFF large 1', 'OFF large 2', 'OFF large 3'}; % visionID
run_opts.cell_type = cell_specification;

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_opts.filepath= ['/Users/colleen/Desktop/Large Cell Summary/', run_opts.date, '/', run_opts.concatname, '/'];


% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames =1:size(datarun.vision.timecourses(1).r,1);% have to input as a vector list of frames, not the number of frames total
datarun=load_data(datarun,opt);
% datarun = load_params(datarun,'verbose',1);
% datarun = load_neurons(datarun);
% datarun = load_sta(datarun);
% datarun = set_polarities(datarun);

if size(cell_specification,2) == 1
    
    [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
    cell_ids=datarun.cell_ids(cell_indices);
    
    
    N = length(cell_ids);
    
    if N> 10
        N = 10;
    end
    
    fig = figure;
    % set(fig, 'PaperSize', [34 7])
    set(fig, 'Position', [0 0 700+200*N 700])
    set(0,'DefaultAxesFontSize',10)
    
    
    % Mosaic
    rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
    
    
    %% Mosaic
    figure(fig)
    subplot(3,2+N, [1,2,2+N+1,2+N+2])
    
    
    rf_c = get(rf, 'Children');
    rf_x = get(rf_c,'Xdata'); % obtain the XData
    rf_y = get(rf_c,'Ydata'); % obtain the YData
    if N>1
        plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
    else
        plot(rf_x, rf_y, 'k')
    end
    
    if strcmp(datarun.stimulus.independent, 'nil')
        colormap gray
    end
    
    axis equal
    axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
    set(gca, 'ydir', 'reverse')
    
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    %%  Timecourses
    params.normalize = 1;
    t = plot_all_timecourses(datarun, cell_specification, params, run_opts);
    
    %% Timecourse
    figure(fig)
    subplot(3,2+N, [2*(2+N)+1,2*(2+N)+2])
    
    t_c = get(t, 'Children');
    xtick = get(t, 'xtick');
    xticklabel = get(t, 'xticklabel');
    
    t_x = get(t_c,'Xdata'); % obtain the XData
    t_y = get(t_c,'Ydata'); % obtain the YData
    
    index= zeros(size(t_x, 1),1);
    index(1:4:end) = 1;
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    plot(x((size(x,1)-run_opts.num_frames):size(x,1),:), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'b')
    
    hold on
    index(1:4:end) = 0;
    index(3:4:end) = 1;
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'r')
    
    index(3:4:end) = 0;
    index(2:4:end) = 1;
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'g')
    
    plot(get(gca, 'xlim'), zeros(1,2), 'k-')
    
    axis tight
    xlabel('time (ms)')
    set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
    set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    %% STAs
    params.padding = 7;
    p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
    
    %% STAs
    figure(fig)
    
    for i = 1:length(cell_ids)
        subplot(3,2+N,2+i)
        
        sta_c = get(p_sta{i}, 'Children');
        sta_image = get(sta_c(3),'Cdata');
        imagesc(sta_image);
        hold on;
        sta_1x = get(sta_c(1),'Xdata');
        sta_1y = get(sta_c(1),'Ydata');
        plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
        
        
        sta_2x = get(sta_c(2),'Xdata');
        sta_2y = get(sta_c(2),'Ydata');
        hold on;
        plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
        axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
        axis square
        axis off
        if i == 10
            break;
        end
        
    end
    
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    % Gen Signal
    gs = computeAndPlot_genSignal(datarun, cell_specification, run_opts );
    
    %% Generator Signal
    figure(fig)
    for i = 1:length(cell_ids)
        subplot(3,2+N, 2+N+2+i)
        
        gs_c = get(gs{i}, 'Children');
        gs_x = get(gs_c,'Xdata'); % obtain the XData
        gs_y = get(gs_c,'Ydata'); % obtain the YData
        plot(gs_x, gs_y, 'o-')
        xlabel('Generator Signal')
        set(gca, 'xlim', [-1 1]);
        if i == 10
            break;
        end
    end
    
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    
    %% ACF
    
    for i = 1:length(cell_ids)
        cell_id = cell_ids(i);
        [time, ccf(i,:)] = plot_ccf(datarun, cell_id);
        subplot(3,2+N, 2+N+2+N +2 +i)
        
        plot(1000*time(time>0), ccf(i, time>0))
        xlabel('time (ms)')
        if i == 10
            break;
        end
    end
    %  avg_ccf = mean(ccf,1);
    %  figure
    % plot(1000*time(time>0), avg_ccf(1, time>0))
    
    
    suptitle({[run_opts.date, run_opts.concatname]; cell_specification{1}})
    if ~exist(run_opts.filepath)
        mkdir(run_opts.filepath)
    end
    export_fig([run_opts.filepath, cell_specification{1}], '-pdf')
    
    
else % More than one cell type, plot them overlayed
    fig = figure;
    % set(fig, 'PaperSize', [34 7])
    set(fig, 'Position', [0 0 700 700])
    set(0,'DefaultAxesFontSize',10)
    Legend= cell(size(cell_specification,2),1);
    
    for i = 1:size(cell_specification,2)
        cell_spec = cell_specification{i};
        
        cmap = hsv(size(cell_specification,2));
        [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_spec);
        cell_ids=datarun.cell_ids(cell_indices);
        N = length(cell_ids);
        %     [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
        % cell_ids=datarun.cell_ids(cell_indices);
        
        
        
        
        
        % Mosaic
        rf = plot_rf_summaries(datarun, cell_spec, 'foa', 0);
        
        
        %% Mosaic
        figure(fig)
        hold on
        
        
        subplot(3,2, [1,2,2+1,2+2])
        
        
        rf_c = get(rf, 'Children');
        rf_x = get(rf_c,'Xdata'); % obtain the XData
        rf_y = get(rf_c,'Ydata'); % obtain the YData
        if N>1
            h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:));
        else
            h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:));
        end
        Legend{i} = cell_spec;
        if i == size(cell_specification,2)
            for j = 1:size(cell_specification,2)
                inter(j) = h{j}(1);
            end
            
            legend_handle = legend([inter] , Legend);
            
        end
        
        hold on
        
        % if strcmp(datarun.stimulus.independent, 'nil')
        %     colormap gray
        % end
        
        axis equal
        axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
        set(gca, 'ydir', 'reverse')
        set(fig, 'HandleVisibility', 'off');
        close all;
        set(fig, 'HandleVisibility', 'on');
        
        
        
        
    end
    
    
    
    
    %% Timecourse
    
    for i = 1:size(cell_specification,2)
        cell_spec = cell_specification{i};
        
        params.normalize = 1;
        t = plot_all_timecourses(datarun, cell_spec, params, run_opts);
        
        %% Timecourse
        figure(fig)
        subplot(3,2, [2*(2)+1,2*(2)+2])
        
        t_c = get(t, 'Children');
        xtick = get(t, 'xtick');
        xticklabel = get(t, 'xticklabel');
        
        t_x = get(t_c,'Xdata'); % obtain the XData
        t_y = get(t_c,'Ydata'); % obtain the YData
        
        index= zeros(size(t_x, 1),1);
        index(2:4:end) = 1;
        x = cell2mat(t_x(logical(index)))';
        y = cell2mat(t_y(logical(index)))';
        H{i} =  plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:),'color', cmap(i,:));
        hold on
        plot(get(gca, 'xlim'), zeros(1,2), 'k-')
        
        axis tight
        xlabel('time (ms)')
        set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
        set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
        set(fig, 'HandleVisibility', 'off');
        close all;
        set(fig, 'HandleVisibility', 'on');
        
        Legend{i} = cell_spec;
        if i == size(cell_specification,2)
            for j = 1:size(cell_specification,2)
                inter(j) = H{j}(1);
            end
            
            legend_handle = legend([inter] , Legend, 'Location', 'Southwest');
            
        end
        
        
    end
    spaces = strfind(cell_specification{2}, ' ');
    
    title_phrase = cell_specification{2}(1:(spaces(2)-1));
    suptitle({[run_opts.date, run_opts.concatname]; title_phrase})
    if ~exist(run_opts.filepath)
        mkdir(run_opts.filepath)
    end
    
    
    export_fig([run_opts.filepath, title_phrase], '-pdf')
end


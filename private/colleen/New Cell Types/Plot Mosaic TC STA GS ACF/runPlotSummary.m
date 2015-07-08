% runPlotSummary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function can create two different types of plots.
% The first one is a plot of one cell type and it's mosaic,
% overlaid timecourses, STAs, ACFs, and nonlinearity for the first 10 or
% fewer cells. The second one is a plot of multiple cell types with the
% mosaics and timecourses (green channel) overlaid. The output of this script
% can be found in /Desktop/Large Cell Summary.

% Functions it calls:
% computeAndPlot_genSignal.m
% plot_all_timecourses.m
% plot_zoomed_STAs.m
% plot_rf_summaries.m
% plot_ccf.m
% export_fig.m
% load_data.m
% get_cell_indices

% Colleen Rhoades
% 2015-07-08
% rhoades@stanford.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_opts.date='2006-06-06-2/'; % one slash at the end
run_opts.concatname='data003-nwpca-cr'; % Name (or modified name) of run, no slashes

% Sometimes the data has two versions of the concate name
% run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];

% Full path to movie xml
run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-4-0.48-33333.xml';


run_opts.save_location_root = '/Users/colleen/Desktop/Large Cell Summary/';
% Number of frames to use for generator signal as well as number of frames
% of the timecourse to display
run_opts.num_frames = 19;

% Number of bins to use for the nonlinearity graph
run_opts.num_bins = 10;

% How much padding to use for zooming in on STAs
params.padding = 7;

% Cell specification can be one cell type or multiple in a cell array.
% Use the same spelling/capitalization as the vision params file
% OFF large-2 in vision = OFF large 2
cell_specification = {'ON large 1','ON large 2','ON large 3','ON large 4','ON large 5'};



%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where to save
run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];

% Used for labeling plot
run_opts.cell_type = cell_specification;

% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
datarun=load_data(datarun,opt);


if size(cell_specification,2) == 1
    %%%%%%%%%%%%%%%%%%%%%%% One Cell Type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % To figure out how many cells are in the cell type
    [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
    cell_ids=datarun.cell_ids(cell_indices);
    N = length(cell_ids);
    
    % Display no more than 10 individual cells for space reasons
    if N> 10
        N = 10;
    end
    
    % Final figure where each component will be plotted
    fig = figure;
    % Size depends on number of cell you are plotting
    set(fig, 'Position', [0 0 700+200*N 700])
    set(0,'DefaultAxesFontSize',10)
    
    %%%%%%%%%%%%%%%%%%%%%%% Mosaic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the axes handle for a figure where all the rf fits are plotted
    rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
    
    
    % Add the mosaic to the output figure
    figure(fig) % make the output figure the current gcf
    subplot(3,2+N, [1,2,2+N+1,2+N+2]) % Top left 2x2 subplot
    
    % get the data plotted in plot_rf_summaries
    rf_c = get(rf, 'Children');
    rf_x = get(rf_c,'Xdata');
    rf_y = get(rf_c,'Ydata');
    
    % type of rf_x and rf_y is only a cell array when there is more than
    % one cell in the type
    if N>1
        plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
    else
        plot(rf_x, rf_y, 'k')
    end
    
    
    %     if strcmp(datarun.stimulus.independent, 'nil')
    %         colormap gray
    %     end
    
    axis equal
    % Set the axis of the mosaic to what the stimulus size was
    axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
    set(gca, 'ydir', 'reverse')
    
    % make the output figure invisible, close the figure that came from
    % plot_rf_summaries, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    %%%%%%%%%%%%%%%%%%%%%%% Overlaid Timecourses %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Turns the axes handle to a plot of the overlaid timecourses
    params.normalize = 1;
    t = plot_all_timecourses(datarun, cell_specification, params, run_opts);
    
    % Make the output figure the gcf
    figure(fig)
    subplot(3,2+N, [2*(2+N)+1,2*(2+N)+2]) % Lower left 2x2 plot
    
    % Extract the lines from the plot from plot_all_timecourses
    t_c = get(t, 'Children');
    xtick = get(t, 'xtick');
    xticklabel = get(t, 'xticklabel');
    
    t_x = get(t_c,'Xdata'); % obtain the XData
    t_y = get(t_c,'Ydata'); % obtain the YData
    
    % Plot all the blue channels together
    index= zeros(size(t_x, 1),1);
    index(1:4:end) = 1; % blue information is 1,4,7 etc
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    % Only plot a certain num_frames
    plot(x((size(x,1)-run_opts.num_frames):size(x,1),:), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'b')
    
    hold on
    % Plot all the red channels together
    
    index(1:4:end) = 0;
    index(3:4:end) = 1; % red information is 3,6,9 etc
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    % Only plot a certain num_frames
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'r')
    
    % Plot all the green channels together
    % Plot green last to it shows in front of the other channels
    index(3:4:end) = 0;
    index(2:4:end) = 1;% green information is 2,5,8 etc
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    % Only plot a certain num_frames
    
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'g')
    
    % Plot a line at 0 to help see the zero crossing
    plot(get(gca, 'xlim'), zeros(1,2), 'k-')
    
    axis tight
    xlabel('time (ms)')
    
    % reduce axis labels to just num_frames
    set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
    set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
    
    % make the output figure invisible, close the figure that came from
    % plot_all_timecourses, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    %%%%%%%%%%%%% Individual STAs for 1-10 cells %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get cell array of axes handles to the plots of the stas zoomed in
    % BUG: if the cell is at the edge of the array the location of the
    % scale bar is wrong and the size looks funky
    p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
    
    % Add individual stas to the output figure
    figure(fig)
    
    for i = 1:length(cell_ids)
        subplot(3,2+N,2+i) % First row subplots 3...12
        
        % Get information about the ith cell's STA
        sta_c = get(p_sta{i}, 'Children');
        sta_image = get(sta_c(3),'Cdata');
        % Plot the image before the fit so it is in the background
        imagesc(sta_image);
        hold on;
        
        % Plot the RF fit
        sta_1x = get(sta_c(1),'Xdata');
        sta_1y = get(sta_c(1),'Ydata');
        plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
        
        % Plot the scale bar (set to 15 pixels)
        sta_2x = get(sta_c(2),'Xdata');
        sta_2y = get(sta_c(2),'Ydata');
        hold on;
        plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
        axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
        axis square
        axis off
        
        % Don't plot more than 10 cells (not enough space)
        if i == 10
            break;
        end
        
    end
    
    % make the output figure invisible, close the figure that came from
    % plot_zoomed_STAs, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    %%%%%%%%%%%%% Nonlinearlity for 1-10 cells %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get a cell array of axes handles for cells 1-10's nonlinearity
    gs = computeAndPlot_genSignal(datarun, cell_specification, run_opts );
    
    % Plot nonlinearity on output figure
    figure(fig)
    for i = 1:length(cell_ids)
        subplot(3,2+N, 2+N+2+i) % Row 2, subplot 3...10
        
        % Get the information from the plot produced by
        % computeAndPlot_genSignal
        gs_c = get(gs{i}, 'Children');
        gs_x = get(gs_c,'Xdata');
        gs_y = get(gs_c,'Ydata');
        plot(gs_x, gs_y, 'o-')
        xlabel('Generator Signal')
        
        % This limit does cut off some GS values, but I want all plots to
        % have the same x scale
        set(gca, 'xlim', [-1 1]);
        
        % Don't plot more than 10 cells, not enough space
        if i == 10
            break;
        end
    end
    
    % make the output figure invisible, close the figure that came from
    % computAndPlot_genSignal, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    
    %%%%%%%%%%%%% Autocorrelation function for 1-10 cells %%%%%%%%%%%%%%%%%
    
    % For each cell (1-10) get the autocorrelation function
    for i = 1:length(cell_ids)
        cell_id = cell_ids(i);
        [time, ccf(i,:)] = plot_ccf(datarun, cell_id);
        
        % Plot in the their row of the output figure
        subplot(3,2+N, 2+N+2+N +2 +i)
        
        % Plot in ms and only for positive time
        plot(1000*time(time>0), ccf(i, time>0))
        xlabel('time (ms)')
        
        % don't plot for more than 10 cells, not enough space
        if i == 10
            break;
        end
    end
    
    % Add an overall title to the output figure
    suptitle({[run_opts.date, run_opts.concatname]; cell_specification{1}})
    
    %%%%%%%%%%%%% Save the output figure %%%%%%%%%%%%%%%%%
    
    if ~exist(run_opts.filepath)
        mkdir(run_opts.filepath)
    end
    export_fig([run_opts.filepath, cell_specification{1}], '-pdf')
    
    
else
    %%%%%%%%%%%%%%%%%%%%%%% More than one Cell Type %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate an output figure that will have the overlaid mosaics and
    % timecourses of the inputted cell types
    fig = figure;
    set(fig, 'Position', [0 0 700 700])
    set(0,'DefaultAxesFontSize',10)
    
    % Will have the text for the legend
    Legend= cell(size(cell_specification,2),1);
    
    for i = 1:size(cell_specification,2)
        % For each cell type
        cell_spec = cell_specification{i};
        
        % color to plot this cell type in
        cmap = distinguishable_colors(size(cell_specification,2));
        
        % Now many cells are in this type
        [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_spec);
        cell_ids=datarun.cell_ids(cell_indices);
        N = length(cell_ids);
        
        %%%%%%%%%%%%%%%%%%%%%%% Mosaic of cell type i %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get axes handle to plot of cell type i's mosaic
        rf = plot_rf_summaries(datarun, cell_spec, 'foa', 0);
        
        % make the output figure the gcf
        figure(fig)
        hold on
        
        % Plot if top 2x2 subplot
        subplot(3,2, [1,2,3,4])
        
        % Get information from the plot of the mosaic
        rf_c = get(rf, 'Children');
        rf_x = get(rf_c,'Xdata'); % obtain the XData
        rf_y = get(rf_c,'Ydata'); % obtain the YData
        
        % rf_x and rf_y are a cell only if more than one cell is in the
        % mosaic
        if N>1
            % Plot each cell type in a different color
            h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:));
        else
            h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:));
        end
        
        % Add the name of the cell type to the legend cell array
        Legend{i} = cell_spec;
        
        % When in the last iteration of the loop
        if i == size(cell_specification,2)
            for j = 1:size(cell_specification,2)
                % Get the line property of the first cell of each cell type
                % so that the colors in the legend are correct
                inter(j) = h{j}(1);
            end
            
            % Specify the data and the string for each cell type
            legend_handle = legend([inter] , Legend);
        end
        
        hold on
        axis equal
        % Make the graph the same size as the stimluus was
        axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
        set(gca, 'ydir', 'reverse')
        
        % make the output figure invisible, close the figure that came from
        % plot_rf_summaries, and then turn the output figure back on
        % Pretends tons of figure from being open at the end of the script
        set(fig, 'HandleVisibility', 'off');
        close all;
        set(fig, 'HandleVisibility', 'on');
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%% Timecourse of cell type i %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:size(cell_specification,2)
        % For each cell type, get the overlaid timecourse
        cell_spec = cell_specification{i};
        
        params.normalize = 1;
        t = plot_all_timecourses(datarun, cell_spec, params, run_opts);
        
        % Make the output figure the gcf
        figure(fig)
        subplot(3,2, [5,6]) % Plot 1x2 subplot
        
        % Get information about the timecourse of cell type i
        t_c = get(t, 'Children');
        xtick = get(t, 'xtick');
        xticklabel = get(t, 'xticklabel');
        
        t_x = get(t_c,'Xdata');
        t_y = get(t_c,'Ydata');
        
        % Only plot the green channel to avoid clutter on the graph
        index= zeros(size(t_x, 1),1);
        index(2:4:end) = 1; %Select only the green channel
        x = cell2mat(t_x(logical(index)))';
        y = cell2mat(t_y(logical(index)))';
        
        % Plot the num_frames of only the green channel with each cell type
        % in a different color
        H{i} =  plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:),'color', cmap(i,:));
        hold on
        
        axis tight
        % Plot a line at zero to see the zero crossing
        plot(get(gca, 'xlim'), zeros(1,2), 'k-')
        
        
        xlabel('time (ms)')
        
        % Set the labels for the right num_frames
        set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
        set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
        
        % make the output figure invisible, close the figure that came from
        % plot_all_timecourses, and then turn the output figure back on
        % Pretends tons of figure from being open at the end of the script
        set(fig, 'HandleVisibility', 'off');
        close all;
        set(fig, 'HandleVisibility', 'on');
        
        % Add the name of the cell type to the cell array of legend entries
        Legend{i} = cell_spec;
        
        % if you are on the last iteration of the loop
        if i == size(cell_specification,2)
            % Select the line axes properites of the first cell of each
            % type
            for j = 1:size(cell_specification,2)
                inter(j) = H{j}(1);
            end
            
            % Set the line properties and string properites of the legend
            legend_handle = legend([inter] , Legend, 'Location', 'Southwest');
            
        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Saving the output figure %%%%%%%%%%%%%%%%%%%%%%
    % Name the figure the overall cell type category (ie ON large when the
    % cell_specification is {'ON large 1' 'ON large 2'}
    spaces = strfind(cell_specification{2}, ' ');
    title_phrase = cell_specification{2}(1:(spaces(2)-1));
    suptitle({[run_opts.date, run_opts.concatname]; title_phrase})
    
    % Make the directory to save in if it doesn't exist
    if ~exist(run_opts.filepath)
        mkdir(run_opts.filepath)
    end
    
    % Save the figure to pdf
    export_fig([run_opts.filepath, title_phrase], '-pdf')
end


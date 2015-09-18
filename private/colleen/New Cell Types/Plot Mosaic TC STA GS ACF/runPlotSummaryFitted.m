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

run_opts.date='2015-08-17-5/'; % one slash at the end
run_opts.concatname='d01-29-norefit'; % Name (or modified name) of run, no slashes
run_opts.concatname_parasol='d01-29-norefit'; % Name (or modified name) of run, no slashes

% Sometimes the data has two versions of the concate name
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data014/data014'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
% 
% run_opts.file_name_parasol = [run_opts.date, '/', run_opts.concatname_parasol, '/',  run_opts.concatname_parasol];
run_opts.file_name_parasol = [run_opts.date, '/', run_opts.concatname_parasol, '/',  'data014/data014'];


% Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';


run_opts.load_location = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname,'/data014'];
run_opts.load_location_parasol = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname_parasol, '/data014'];

run_opts.save_location_root = '/Users/colleen/Desktop/Large Cell Summary Fitted/';
% Number of frames to use for generator signal as well as number of frames
% of the timecourse to display
run_opts.num_frames = 20;

% Number of bins to use for the nonlinearity graph
run_opts.num_bins = 10;

% How much padding to use for zooming in on STAs
params.padding = 7;

% Cell specification can be one cell type or multiple in a cell array.
% Use the same spelling/capitalization as the vision params file
% OFF large-2 in vision = OFF large 2
cell_specification = {'ON large 3'};
cell_specification_parasol = {'ON parasol'};
run_parasol = 1;
run_genSignal = 0;
num_electrodes = 512;

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where to save
run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/data014/'];
run_opts.filepath_parasol= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname_parasol, '/data014/'];

% Used for labeling plot
run_opts.cell_type = cell_specification;

% load fitted datadata014
s=dir(fullfile(run_opts.load_location,'*.mat'));

if size(cell_specification,2) > 1
    for cell_spec = 1:size(cell_specification,2)
        cells = cell_specification{cell_spec};
        count = 1;
        for i = 1:size(s,1)
            period = strfind(s(i).name, '.');
            if sum(strcmp(s(i).name(1:period-1), cells))>0
                fitted_cells{cell_spec}(count) =  load([run_opts.load_location, '/', s(i).name]);
                count = count+1;
            end
            
            
        end
    end
    
else
    
    count = 1;
    for i = 1:size(s,1)
        period = strfind(s(i).name, '.');
        if sum(strcmp(s(i).name(1:period-1), cell_specification))>0
            fitted_cells(count) =  load([run_opts.load_location, '/', s(i).name]);
            count = count+1;
        end
        
        
    end
    if run_parasol
        q=dir(fullfile(run_opts.load_location_parasol,'*.mat'));
        count = 1;
        for i = 1:size(q,1)
            period = strfind(q(i).name, '.');
            if sum(strcmp(q(i).name(1:period-1), cell_specification_parasol))>0
                fitted_cells_parasol(count) =  load([run_opts.load_location_parasol, '/', q(i).name]);
                count = count+1;
            end
            
            
        end
    end
end


% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1);
opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
datarun=load_data(datarun,opt);
datarun.names.rrs_ei_path = ['/Volumes/Analysis/', run_opts.file_name, '.ei'];

datarun = load_ei(datarun, cell_specification, 'array_type', num_electrodes);


% Load Data
datarun_parasol.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name_parasol, '.neurons'];
datarun_parasol.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name_parasol, '.params'];
datarun_parasol.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name_parasol, '.sta'];
% datarun.names.rrs_ei_path = ['/Volumes/Analysis/', run_opts.file_name_parasol, '.ei'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
datarun_parasol=load_data(datarun_parasol,opt);
% datarun_parasol = load_ei(datarun_parasol, cell_specification_parasol);



if size(cell_specification,2) == 1
    %%%%%%%%%%%%%%%%%%%%%%% One Cell Type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % To figure out how many cells are in the cell type
    [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
    cell_ids=datarun.cell_ids(cell_indices);
    %     N = length(cell_ids);
    
    N = size(fitted_cells.prop.sta_axes,2);
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
    
    
    % Add the mosaic to the output figure
    figure(fig) % make the output figure the current gcf
    subplot(3,2+N, [1,2]) % Top left 2x2 subplot
    
    % Get the axes handle for a figure where all the rf fits are plotted
    rf_c = get(fitted_cells.prop.mosaic_axes, 'Children');
    rf_x = get(rf_c,'Xdata');
    rf_y = get(rf_c,'Ydata');
    
    
    
    
    
    
    % type of rf_x and rf_y is only a cell array when there is more than
    % one cell in the type
    if N>1
        plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
    else
        plot(rf_x, rf_y, 'k')
    end
    
    axis equal
    % Set the axis of the mosaic to what the stimulus size was
    %     axis([bounds])
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
    %     params.normalize = 1;
    %     t = plot_all_timecourses(datarun, cell_specification, params, run_opts);
    
    % Make the output figure the gcf
    figure(fig)
    subplot(3,2+N, [(2+N)+1,(2+N)+2]) % Lower left 2x2 plot
    
    % Extract the lines from the plot from plot_all_timecourses
    t_c = get(fitted_cells.prop.t_axes, 'Children');
    
    
    
    %     xtick = get(t, 'xtick');
    %     xticklabel = get(t, 'xticklabel');
    
    t_x = get(t_c,'Xdata'); % obtain the XData
    t_y = get(t_c,'Ydata'); % obtain the YData
    time = round(fliplr(0:-datarun.stimulus.interval*8.33:-(run_opts.num_frames-1)*datarun.stimulus.interval*8.33));
    
    for j = 1:size(t_x,1)
        try
            plot(time, t_y{j}(end-run_opts.num_frames+1:end), 'color', get(t_c(j), 'color'))
            hold on
        catch
            plot(time, t_y(end-run_opts.num_frames+1:end), 'color', get(t_c(j), 'color'))
            hold on
        end
        
    end
    
    plot(time, fitted_cells.prop.avg_timecourse(end-run_opts.num_frames+1:end), '-k')
    set(gca, 'xtick', time(1:2:end))
    set(gca, 'xticklabel', time(1:2:end));
    
    
    % Plot a line at 0 to help see the zero crossing
    plot(get(gca, 'xlim'), zeros(1,2), 'k-')
    xlabel('time (ms)')
    
    %     axis tight
    %     xlabel('time (ms)')
    %
    %     % reduce axis labels to just num_frames
    %     set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
    %     set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
    %
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
    %     p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
    %% new
    stas = fitted_cells.prop.sta_axes;
    
    %%
    % Add individual stas to the output figure
    figure(fig)
    
    for i = 1:length(stas)
        subplot(3,2+N,2+i) % First row subplots 3...12
        
        
        % Get information about the ith cell's STA
        
        sta_c = get(stas{i}, 'Children');
        sta_image = get(sta_c(2),'cdata');
        % Plot the image before the fit so it is in the background
        imagesc(sta_image);
        hold on;
        
        % Plot the RF fit
        sta_1x = get(sta_c(1),'Xdata');
        sta_1y = get(sta_c(1),'Ydata');
        
        
        
        m_x =  mean(sta_1x);
        m_y =  mean(sta_1y);
        sd_x =  std(sta_1x);
        sd_y =  std(sta_1y);
        pad_factor = params.padding;
        scale = [1 1];
        aspect_ratio = 1;
        
        xdiff = sd_x / 2 * pad_factor;
        ydiff = sd_y / 2 * pad_factor;
        ydiff = max(ydiff, xdiff/aspect_ratio);
        xdiff = max(xdiff, ydiff*aspect_ratio);
        
        xstart = m_x - xdiff;
        xend   = m_x + xdiff;
        ystart = m_y - ydiff;
        yend   = m_y + ydiff;
        bounds = [scale(1) scale(1) scale(2) scale(2)].*[xstart xend ystart yend];
        
        
        
        
        plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
        
        % Plot the scale bar (set to 15 pixels)
        %         sta_2x = get(sta_c(2),'Xdata');
        %         sta_2y = get(sta_c(2),'Ydata');
        %         hold on;
        %         plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
        %         axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
        axis([bounds])
        axis square
        axis off
        %         axis equal
        % Set the axis of the mosaic to what the stimulus size was
        %     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
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
    
    if run_genSignal ==1
        % Get a cell array of axes handles for cells 1-10's nonlinearity
        [gs, slope] = computeAndPlot_genSignal(datarun, cell_specification, run_opts );
        if run_parasol
            [gs_p, slope_p] = computeAndPlot_genSignal(datarun_parasol, cell_specification_parasol, run_opts );
        end
    end
    
    %% EI
    figure(fig);
    array_info = load_array_info(datarun,2);
    position= tformfwd(array_info.T_array_to_vision_ei,datarun.ei.position);
    %     positions(:,1)=datarun.stimulus.field_width-positions(:,1)+.5;
    %     position(:,2)=datarun.stimulus.field_height-position(:,2);
    elec_handles = cell(min(length(cell_ids),10),1);
    for i = 1:length(cell_ids)
        subplot(3,2+N, 2+N+2+i) % Row 2, subplot 3...10
        
        elec_handles{i} = plot_ei(datarun, cell_ids(i));
        
        
        for j = 1:size(elec_handles,2)
            try
                fac = get(elec_handles{i}(j), 'faces');
                vert = get(elec_handles{i}(j), 'vertices');
                patch('faces', fac, 'vertices', vert)
                hold on
            catch
            end
        end
        set(gca, 'ydir','reverse', 'xdir', 'reverse')
        axis equal;
        axis tight;
        box on;
        set(gca, 'xtick', [], 'ytick', [])
        %         axis off;
        
        range=[min(position(:,1)) max(position(:,1)) min(position(:,2)) max(position(:,2))*1.1];
        axis(range);
        if i == 10
            break;
        end
    end
    
    
    %     figure;
    %     plot(elec_handles)
    % Plot nonlinearity on output figure
    %replace with EIs
    %     figure(fig)
    %     for i = 1:length(cell_ids)
    %         subplot(3,2+N, 2+N+2+i) % Row 2, subplot 3...10
    %
    %         % Get the information from the plot produced by
    %         % computeAndPlot_genSignal
    %         gs_c = get(gs{i}, 'Children');
    %         gs_x = get(gs_c,'Xdata');
    %         gs_y = get(gs_c,'Ydata');
    %         plot(gs_x, gs_y, 'o-')
    %         xlabel('Generator Signal')
    %
    %         % This limit does cut off some GS values, but I want all plots to
    %         % have the same x scale
    %         set(gca, 'xlim', [-1 1]);
    %
    %         % Don't plot more than 10 cells, not enough space
    %         if i == 10
    %             break;
    %         end
    %     end
    
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
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RF = [num2str(mean(fitted_cells.prop.RF_diameter)), ' +/- ', num2str(std(fitted_cells.prop.RF_diameter))];
    FR = [num2str(mean(fitted_cells.prop.FR)), ' +/- ', num2str(std(fitted_cells.prop.FR)), ' Hz'];
    amp= num2str(fitted_cells.prop.amplitude);
    zc = num2str(fitted_cells.prop.zc_time);
    peak = num2str(fitted_cells.prop.peak_time);
    trough = num2str(fitted_cells.prop.trough_time);
    if run_genSignal == 1
        slope_str = [num2str(mean(slope)), ' +/- ', num2str(std(slope))];
        slope_str_p = [num2str(mean(slope_p)), ' +/- ', num2str(std(slope_p))];
        
    end
    
    RF_p = [num2str(mean(fitted_cells_parasol.prop.RF_diameter)), ' +/- ', num2str(std(fitted_cells_parasol.prop.RF_diameter)), ' ', char(0181), 'm'];
    FR_p = [num2str(mean(fitted_cells_parasol.prop.FR)), ' +/- ', num2str(std(fitted_cells_parasol.prop.FR)), ' Hz'];
    amp_p= num2str(fitted_cells_parasol.prop.amplitude);
    zc_p = num2str(fitted_cells_parasol.prop.zc_time);
    peak_p = num2str(fitted_cells_parasol.prop.peak_time);
    trough_p = num2str(fitted_cells_parasol.prop.trough_time);
    
    figure(fig)
    % create the data
    % Create the column and row names in cell arrays
    % cnames = {'1','2','3','4','5','6','7','8'};
    if run_genSignal == 1
        rnames = {'RF diameter','FR', 'STA peak amp','Time to Zero','Time to Peak','Time to Trough','Slope of SNL at 0'};
    else
        rnames = {'RF diameter','FR', 'STA peak amp','Time to Zero','Time to Peak','Time to Trough'};
        
    end
    % % Create the uitable
    % % t = uitable(fig,'Data',{rnames{1}, RF, RF_p; rnames{2},FR, FR_p; rnames{3},amp, amp_p; rnames{4},zc, zc_p; rnames{5},peak, peak_p; rnames{6},trough, trough_p; rnames{7},slope_str, slope_str_p},...
    % %             'ColumnName',{[], 'Large', 'Parasol'},...
    % %             'RowName',[],...
    % %             'ColumnWidth',{90, 130, 130});
    %
    %
    %         t = uitable('Data',rand(3),...
    %             'ColumnName',{[], 'Large', 'Parasol'},...
    %             'RowName',[],...
    %             'ColumnWidth',{90, 100, 100}, 'Parent', fig);
    % subplot(3,2+N,[2*(2+N)+1,2*(2+N)+2]),plot(3)
    % pos = get(subplot(3,2+N,[2*(2+N)+1,2*(2+N)+2]),'position');
    % delete(subplot(3,2+N,[2*(2+N)+1,2*(2+N)+2]))
    % set(t,'units','normalized')
    % set(t,'position',pos)
    % axis off;
    
    k = subplot(3,2+N,[2*(2+N)+1,2*(2+N)+2]);%,plot(3)
    % k = subplot(3,2+N,[1,1])%,plot(3)
    
    xl = xlim(k);
    xPos = xl(1) + diff(xl) / 3;
    yl = ylim(k);
    yPos = yl(1) + diff(yl) / 3;
    % first column
    if N > 5
        if run_genSignal == 1
            t1 = text(-0.6, yPos, sprintf('\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', rnames{1}, rnames{2}, rnames{3}, rnames{4}, rnames{5}, rnames{6}, rnames{7}), 'Parent', k, 'fontsize',11);
            t2 = text(xPos-.5, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Large', RF, FR, amp, zc, peak, trough, slope_str), 'Parent', k, 'fontsize',11);
            t3 = text(2*xPos-0.2, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Parasol',RF_p, FR_p, amp_p, zc_p, peak_p, trough_p, slope_str_p), 'Parent', k, 'fontsize',11);
        else
            t1 = text(-0.6, yPos, sprintf('\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%', rnames{1}, rnames{2}, rnames{3}, rnames{4}, rnames{5}, rnames{6}), 'Parent', k, 'fontsize',11);
            t2 = text(xPos-.5, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%', 'Large', RF, FR, amp, zc, peak, trough), 'Parent', k, 'fontsize',11);
            t3 = text(2*xPos-0.2, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Parasol',RF_p, FR_p, amp_p, zc_p, peak_p, trough_p), 'Parent', k, 'fontsize',11);
            
        end
        
    else
        if run_genSignal ==1
            t1 = text(0, yPos, sprintf('\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', rnames{1}, rnames{2}, rnames{3}, rnames{4}, rnames{5}, rnames{6}, rnames{7}), 'Parent', k, 'fontsize',11);
            t2 = text(xPos, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Large', RF, FR, amp, zc, peak, trough, slope_str), 'Parent', k, 'fontsize',11);
            t3 = text(2*xPos, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Parasol',RF_p, FR_p, amp_p, zc_p, peak_p, trough_p, slope_str_p), 'Parent', k, 'fontsize',11);
        else
            
            t1 = text(0, yPos, sprintf('\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', rnames{1}, rnames{2}, rnames{3}, rnames{4}, rnames{5}, rnames{6}), 'Parent', k, 'fontsize',11);
            t2 = text(xPos, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Large', RF, FR, amp, zc, peak, trough), 'Parent', k, 'fontsize',11);
            t3 = text(2*xPos, yPos, sprintf('%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n%s\n\n', 'Parasol',RF_p, FR_p, amp_p, zc_p, peak_p, trough_p), 'Parent', k, 'fontsize',11);
            
        end
        
        
    end
    % % t = uitable(fig,'Data',{rnames{1}, RF, RF_p; rnames{2},FR, FR_p; rnames{3},amp, amp_p; rnames{4},zc, zc_p; rnames{5},peak, peak_p; rnames{6},trough, trough_p; rnames{7},slope_str, slope_str_p},...
    axis off
    % set(t1,'RightAlignment','center')
    
    % Add an overall title to the output figure
    suptitle({[run_opts.date, run_opts.concatname]; cell_specification{1}})
    
    %%%%%%%%%%%%% Save the output figure %%%%%%%%%%%%%%%%%
    
    if ~exist(run_opts.filepath)
        mkdir(run_opts.filepath)
    end
    export_fig([run_opts.filepath, cell_specification{1}], '-pdf')
    
    %% More than one cell type
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Add the mosaic to the output figure
        figure(fig) % make the output figure the current gcf
        hold on;
        subplot(3,2, [1,2,3,4])
        % Get the axes handle for a figure where all the rf fits are plotted
        rf_c = get(fitted_cells{i}.prop.mosaic_axes, 'Children');
        rf_x = get(rf_c,'Xdata');
        rf_y = get(rf_c,'Ydata');
        
        
        
        
        
        
        % type of rf_x and rf_y is only a cell array when there is more than
        % one cell in the type
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
        % Set the axis of the mosaic to what the stimulus size was
        %     axis([bounds])
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
        %         t = plot_all_timecourses(datarun, cell_spec, params, run_opts);
        % Make the output figure the gcf
        figure(fig)
        subplot(3,2, [5,6]) % Plot 1x2 subplot
        
        t_c = get(fitted_cells{i}.prop.t_axes, 'Children');
        
        %         t_c = get(t, 'Children');
        %         xtick = get(t, 'xtick');
        %         xticklabel = get(t, 'xticklabel');
        
        t_x = get(t_c,'Xdata');
        t_y = get(t_c,'Ydata');
        time = round(fliplr(0:-datarun.stimulus.interval*8.33:-(run_opts.num_frames-1)*datarun.stimulus.interval*8.33));
        for j = 1:size(t_x,1)
            
            if get(t_c(j), 'color') == [0,0,0]
                try
                    H{i} = plot(time, t_y{j}(end-run_opts.num_frames+1:end), 'color', cmap(i,:));
                    hold on
                catch
                    H{i} =plot(time, t_y(end-run_opts.num_frames+1:end), 'color',cmap(i,:));
                    hold on
                end
                
            end
            
        end
        
        %         H{i} =  plot(time, fitted_cells.prop.avg_timecourse(end-run_opts.num_frames+1:end), 'k')
        set(gca, 'xtick', time(1:2:end))
        set(gca, 'xticklabel', time(1:2:end));
        
        axis tight
        
        hold on
        
        
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
    plot(get(gca, 'xlim'), zeros(1,2), 'k-')
    xlabel('time (ms)')
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
    
    
    
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
    export_fig([run_opts.filepath, [title_phrase, ' Summary']], '-pdf')
end


% 
% run_opts.date='2009-12-03-0/'; % one slash at the end
% run_opts.concatname='data003-nwpca'; % Name (or modified name) of run, no slashes
% 
% % Sometimes the data has two versions of the concate name
% % run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'daat003/daat003'];
% 
% % Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111-64x48.xml';
% 
% 
% run_opts.save_location_root = '/Users/colleen/Desktop/Berstein Talk/Mosaics/';
% % Number of frames to use for generator signal as well as number of frames
% % of the timecourse to display
% run_opts.num_frames = 19;
% 
% % Number of bins to use for the nonlinearity graph
% run_opts.num_bins = 10;
% 
% % How much padding to use for zooming in on STAs
% params.padding = 7;
% 
% % Cell specification can be one cell type or multiple in a cell array.
% % Use the same spelling/capitalization as the vision params file
% % OFF large-2 in vision = OFF large 2
% cell_specification = {'ON parasol'};
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where to save
% run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];
% 
% % Used for labeling plot
% run_opts.cell_type = cell_specification;
% 
% % Load Data
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
% opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
% datarun=load_data(datarun,opt);
% 
% filepath=[sprintf('%s%s%s.pdf',[run_opts.save_location_root,'/',run_opts.date,'/',run_opts.concatname, '/'])];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end
% 
% 
%     % To figure out how many cells are in the cell type
%     [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
%     cell_ids=datarun.cell_ids(cell_indices);
%     N = length(cell_ids);
%     
%     % Display no more than 10 individual cells for space reasons
%     if N> 5
%         N = 5;
%     end
%     
%     % Final figure where each component will be plotted
%     fig = figure;
%     % Size depends on number of cell you are plotting
%     set(fig, 'Position', [0 0 700 700])
%     set(0,'DefaultAxesFontSize',10)
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Mosaic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Get the axes handle for a figure where all the rf fits are plotted
%     rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
%     
%         % Add the mosaic to the output figure
%     figure(fig) % make the output figure the current gcf
%     subplot(3,2+N, [1,2,2+N+1,2+N+2]) % Top left 2x2 subplot
%     
%     % get the data plotted in plot_rf_summaries
%     rf_c = get(rf, 'Children');
%     rf_x = get(rf_c,'Xdata');
%     rf_y = get(rf_c,'Ydata');
%     
%     % type of rf_x and rf_y is only a cell array when there is more than
%     % one cell in the type
%     if N>1
%         plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
%     else
%         plot(rf_x, rf_y, 'k')
%     end
%     
%     
%     %     if strcmp(datarun.stimulus.independent, 'nil')
%     %         colormap gray
%     %     end
%     
%     axis equal
%     % Set the axis of the mosaic to what the stimulus size was
%     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
%     set(gca, 'ydir', 'reverse')
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_rf_summaries, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
% %     set(fig, 'HandleVisibility', 'off');
% %     close all;
% %     set(fig, 'HandleVisibility', 'on');
% 
% 
%  % Get cell array of axes handles to the plots of the stas zoomed in
%     % BUG: if the cell is at the edge of the array the location of the
%     % scale bar is wrong and the size looks funky
%     p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
%     
%     % Add individual stas to the output figure
%     figure(fig)
%     
%     for i = 1:length(cell_ids)
%         subplot(3,2+N,2+i) % First row subplots 3...12
%         
%         % Get information about the ith cell's STA
%         sta_c = get(p_sta{i}, 'Children');
%         sta_image = get(sta_c(3),'Cdata');
%         % Plot the image before the fit so it is in the background
%         imagesc(sta_image);
%         hold on;
%         
%         % Plot the RF fit
%         sta_1x = get(sta_c(1),'Xdata');
%         sta_1y = get(sta_c(1),'Ydata');
%         plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
%         
%         % Plot the scale bar (set to 15 pixels)
%         sta_2x = get(sta_c(2),'Xdata');
%         sta_2y = get(sta_c(2),'Ydata');
%         hold on;
%         plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
%         axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
%         axis square
%         axis off
%         
%         % Don't plot more than 10 cells (not enough space)
%         if i == 5
%             break;
%         end
%         
%     end
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_zoomed_STAs, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
% 
%     hgexport(fig, [filepath,cell_specification{1}])

%     
%     
%     %% New cell type
%     
%     clear
% run_opts.date='2005-08-03-0/'; % one slash at the end
% run_opts.concatname='data001'; % Name (or modified name) of run, no slashes
% 
% % Sometimes the data has two versions of the concate name
% % run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data001'];
% 
% % Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
% 
% 
% run_opts.save_location_root = '/Users/colleen/Desktop/Berstein Talk/Mosaics/';
% % Number of frames to use for generator signal as well as number of frames
% % of the timecourse to display
% run_opts.num_frames = 19;
% 
% % Number of bins to use for the nonlinearity graph
% run_opts.num_bins = 10;
% 
% % How much padding to use for zooming in on STAs
% params.padding = 10;
% 
% % Cell specification can be one cell type or multiple in a cell array.
% % Use the same spelling/capitalization as the vision params file
% % OFF large-2 in vision = OFF large 2
% cell_specification = {'ON midget'};
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where to save
% run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];
% 
% % Used for labeling plot
% run_opts.cell_type = cell_specification;
% 
% % Load Data
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
% opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
% datarun=load_data(datarun,opt);
% 
% filepath=[sprintf('%s%s%s.pdf',[run_opts.save_location_root,'/',run_opts.date,'/',run_opts.concatname, '/'])];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end
% 
% 
%     % To figure out how many cells are in the cell type
%     [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
%     cell_ids=datarun.cell_ids(cell_indices);
%     N = length(cell_ids);
%     
%     % Display no more than 10 individual cells for space reasons
%     if N> 5
%         N = 5;
%     end
%     
%     % Final figure where each component will be plotted
%     fig = figure;
%     % Size depends on number of cell you are plotting
%     set(fig, 'Position', [0 0 700 700])
%     set(0,'DefaultAxesFontSize',10)
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Mosaic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Get the axes handle for a figure where all the rf fits are plotted
%     rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
%     
%         % Add the mosaic to the output figure
%     figure(fig) % make the output figure the current gcf
%     subplot(3,2+N, [1,2,2+N+1,2+N+2]) % Top left 2x2 subplot
%     
%     % get the data plotted in plot_rf_summaries
%     rf_c = get(rf, 'Children');
%     rf_x = get(rf_c,'Xdata');
%     rf_y = get(rf_c,'Ydata');
%     
%     % type of rf_x and rf_y is only a cell array when there is more than
%     % one cell in the type
%     if N>1
%         plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
%     else
%         plot(rf_x, rf_y, 'k')
%     end
%     
%     
%     %     if strcmp(datarun.stimulus.independent, 'nil')
%     %         colormap gray
%     %     end
%     
%     axis equal
%     % Set the axis of the mosaic to what the stimulus size was
%     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
%     set(gca, 'ydir', 'reverse')
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_rf_summaries, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
% %     set(fig, 'HandleVisibility', 'off');
% %     close all;
% %     set(fig, 'HandleVisibility', 'on');
% 
% 
%  % Get cell array of axes handles to the plots of the stas zoomed in
%     % BUG: if the cell is at the edge of the array the location of the
%     % scale bar is wrong and the size looks funky
%     p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
%     
%     % Add individual stas to the output figure
%     figure(fig)
%     
%     for i = 10:length(cell_ids)
%         subplot(3,2+N,2+i) % First row subplots 3...12
%         
%         % Get information about the ith cell's STA
%         sta_c = get(p_sta{i}, 'Children');
%         sta_image = get(sta_c(3),'Cdata');
%         % Plot the image before the fit so it is in the background
%         imagesc(sta_image);
%         hold on;
%         
%         % Plot the RF fit
%         sta_1x = get(sta_c(1),'Xdata');
%         sta_1y = get(sta_c(1),'Ydata');
%         plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
%         
%         % Plot the scale bar (set to 15 pixels)
%         sta_2x = get(sta_c(2),'Xdata');
%         sta_2y = get(sta_c(2),'Ydata');
%         hold on;
%         plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
%         axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
%         axis square
%         axis off
%         
%         % Don't plot more than 10 cells (not enough space)
%         if i == 5
%             break;
%         end
%         
%     end
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_zoomed_STAs, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
% 
%     hgexport(fig, [filepath,cell_specification{1}])

    %% Five cell types
%     clear
%     run_opts.date='2013-05-28-4/'; % one slash at the end
% run_opts.concatname='data000'; % Name (or modified name) of run, no slashes
% 
% % Sometimes the data has two versions of the concate name
% % run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data000'];
% 
% % Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
% 
% 
% run_opts.save_location_root = '/Users/colleen/Desktop/Berstein Talk/Mosaics/';
% % Number of frames to use for generator signal as well as number of frames
% % of the timecourse to display
% run_opts.num_frames = 19;
% 
% % Number of bins to use for the nonlinearity graph
% run_opts.num_bins = 10;
% 
% % How much padding to use for zooming in on STAs
% params.padding = 7;
% 
% % Cell specification can be one cell type or multiple in a cell array.
% % Use the same spelling/capitalization as the vision params file
% % OFF large-2 in vision = OFF large 2
% cell_specification = {'ON large 1'};
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where to save
% run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];
% 
% % Used for labeling plot
% run_opts.cell_type = cell_specification;
% 
% % Load Data
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
% opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
% datarun=load_data(datarun,opt);
% 
% filepath=[sprintf('%s%s%s.pdf',[run_opts.save_location_root,'/',run_opts.date,'/',run_opts.concatname, '/'])];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end
% 
% 
%     % To figure out how many cells are in the cell type
%     [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
%     cell_ids=datarun.cell_ids(cell_indices);
%     N = length(cell_ids);
%     
%     % Display no more than 10 individual cells for space reasons
%     if N> 5
%         N = 5;
%     end
%     
%     % Final figure where each component will be plotted
%     fig = figure;
%     % Size depends on number of cell you are plotting
%     set(fig, 'Position', [0 0 1000 300])
%     set(0,'DefaultAxesFontSize',10)
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Mosaic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Get the axes handle for a figure where all the rf fits are plotted
%     rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
%     
%         % Add the mosaic to the output figure
%     figure(fig) % make the output figure the current gcf
%     subplot(3,2+N, [1,2,2+N+1,2+N+2]) % Top left 2x2 subplot
%     
%     % get the data plotted in plot_rf_summaries
%     rf_c = get(rf, 'Children');
%     rf_x = get(rf_c,'Xdata');
%     rf_y = get(rf_c,'Ydata');
%     
%     % type of rf_x and rf_y is only a cell array when there is more than
%     % one cell in the type
%     if N>1
%         plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
%     else
%         plot(rf_x, rf_y, 'k')
%     end
%     
%     
%     %     if strcmp(datarun.stimulus.independent, 'nil')
%     %         colormap gray
%     %     end
%     
%     axis equal
%     % Set the axis of the mosaic to what the stimulus size was
% %     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
%     axis([ 20 50 0 25])
% 
%     set(gca, 'ydir', 'reverse')
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_rf_summaries, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
% %     set(fig, 'HandleVisibility', 'off');
% %     close all;
% %     set(fig, 'HandleVisibility', 'on');
% 
% 
%  %%%%%%%%%%%%%%%%%%%%%%% Overlaid Timecourses %%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Turns the axes handle to a plot of the overlaid timecourses
%     params.normalize = 1;
%     t = plot_all_timecourses(datarun, cell_specification, params, run_opts);
%     
%     % Make the output figure the gcf
%     figure(fig)
%     subplot(3,2+N, [2*(2+N)+1,2*(2+N)+2]) % Lower left 2x2 plot
%     
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     xtick = get(t, 'xtick');
%     xticklabel = get(t, 'xticklabel');
%     
%     t_x = get(t_c,'Xdata'); % obtain the XData
%     t_y = get(t_c,'Ydata'); % obtain the YData
%     
%     % Plot all the blue channels together
%     index= zeros(size(t_x, 1),1);
%     index(1:4:end) = 1; % blue information is 1,4,7 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1),:), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [0 0 1 0.05])
%     
%     hold on
%     % Plot all the red channels together
%     
%     index(1:4:end) = 0;
%     index(3:4:end) = 1; % red information is 3,6,9 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [1 0 0 0.05])
%     
%     % Plot all the green channels together
%     % Plot green last to it shows in front of the other channels
%     index(3:4:end) = 0;
%     index(2:4:end) = 1;% green information is 2,5,8 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:),'Color', [0 1 0 0.05])
%     
%     % Plot a line at 0 to help see the zero crossing
%     plot(get(gca, 'xlim'), zeros(1,2), 'k-')
%     
%     axis tight
%     xlabel('time (ms)')
%     
%     % reduce axis labels to just num_frames
%     set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
%     set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_all_timecourses, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
%     
%     %%%%%% Individual STAs%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  % Get cell array of axes handles to the plots of the stas zoomed in
%     % BUG: if the cell is at the edge of the array the location of the
%     % scale bar is wrong and the size looks funky
%     p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
%     
%     % Add individual stas to the output figure
%     figure(fig)
%     
%     for i = 1:length(cell_ids)
%         subplot(3,2+N,2+i) % First row subplots 3...12
%         
%         % Get information about the ith cell's STA
%         sta_c = get(p_sta{i}, 'Children');
%         sta_image = get(sta_c(3),'Cdata');
%         % Plot the image before the fit so it is in the background
%         imagesc(sta_image);
%         hold on;
%         
%         % Plot the RF fit
%         sta_1x = get(sta_c(1),'Xdata');
%         sta_1y = get(sta_c(1),'Ydata');
%         plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
%         
%         % Plot the scale bar (set to 15 pixels)
%         sta_2x = get(sta_c(2),'Xdata');
%         sta_2y = get(sta_c(2),'Ydata');
%         hold on;
%         plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
%         axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
%         axis square
%         axis off
%         
%         % Don't plot more than 10 cells (not enough space)
%         if i == 5
%             break;
%         end
%         
%     end
%     
%         % make the output figure invisible, close the figure that came from
%     % plot_zoomed_STAs, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
%     
%     %%%%%%%%%%%%%%%% PCs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     [cell_indicesSBC] = get_cell_indices(datarun, {'SBC'});
%     cell_idsSBC=datarun.cell_ids(cell_indicesSBC);
%     
%       [cell_indicesONPar] = get_cell_indices(datarun, {'ON parasol'});
%     cell_idsONPar=datarun.cell_ids(cell_indicesONPar);
%     
%       [cell_indicesOFFPar] = get_cell_indices(datarun, {'OFF parasol'});
%     cell_idsOFFPar=datarun.cell_ids(cell_indicesOFFPar);
%     
%       [cell_indicesONMid] = get_cell_indices(datarun, {'ON midget'});
%     cell_idsONMid=datarun.cell_ids(cell_indicesONMid);
%     
%       [cell_indicesOFFMid] = get_cell_indices(datarun, {'OFF midget'});
%     cell_idsOFFMid=datarun.cell_ids(cell_indicesOFFMid);
%    
%     
%           [cell_indicesONLar] = get_cell_indices(datarun, {'ON large 1'});
%     cell_idsONLar=datarun.cell_ids(cell_indicesONLar);
%     
%     N = length(cell_ids);
%     
%        areaSBC = zeros(length(cell_indicesSBC),1);
%        areaONPar = zeros(length(cell_indicesONPar),1);
%        areaOFFPar = zeros(length(cell_indicesOFFPar),1);
%        areaONMid = zeros(length(cell_indicesONMid),1);
%        areaOFFMid = zeros(length(cell_indicesOFFMid),1);
%        areaONLar = zeros(length(cell_indicesONLar),1);
% 
% %        numON = length(areaONPar) + length(areaONMid);
% %        numOFF = length(areaOFFPar) + length(areaOFFMid);
% %        
% %        scale_OFF_by = numON/numOFF;
%        for i = 1:length(cell_indicesSBC)
%            fit = datarun.stas.fits{cell_indicesSBC(i)};
%            areaSBC(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end
%        
%        for i = 1:length(cell_indicesONPar)
%            fit = datarun.stas.fits{cell_indicesONPar(i)};
%            areaONPar(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end
%        
%        for i = 1:length(cell_indicesOFFPar)
%            fit = datarun.stas.fits{cell_indicesOFFPar(i)};
%            areaOFFPar(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end
%        
%        for i = 1:length(cell_indicesONMid)
%            fit = datarun.stas.fits{cell_indicesONMid(i)};
%            areaONMid(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end
%        
%        for i = 1:length(cell_indicesOFFMid)
%            fit = datarun.stas.fits{cell_indicesOFFMid(i)};
%            areaOFFMid(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end
%     for i = 1:length(cell_indicesONLar)
%            fit = datarun.stas.fits{cell_indicesONLar(i)};
%            areaONLar(i) = pi*fit.sd(1).*fit.sd(2)*datarun.stimulus.stixel_width*5.5; % 5.5 um/pixel;
%        end       
%    
%     
%     t = plot_all_timecourses(datarun, {'SBC'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_ySBC = get(t_c,'Ydata'); % obtain the YData
%  
%         t = plot_all_timecourses(datarun, {'ON parasol'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_yONPar = get(t_c,'Ydata'); % obtain the YData
%     
%         t = plot_all_timecourses(datarun, {'OFF parasol'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_yOFFPar = get(t_c,'Ydata'); % obtain the YData
%     
%         t = plot_all_timecourses(datarun, {'ON midget'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_yONMid = get(t_c,'Ydata'); % obtain the YData
%     
%         t = plot_all_timecourses(datarun, {'OFF midget'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_yOFFMid = get(t_c,'Ydata'); % obtain the YData
%     
%     
%             t = plot_all_timecourses(datarun, {'ON large 1'}, params, run_opts);
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     t_yONLar = get(t_c,'Ydata'); % obtain the YData
%     
%     
%     counter= 1;
%     for i =1:4:size(t_ySBC,1)-2
%         tcSBC(counter,:) = mean([t_ySBC{i}; t_ySBC{i+1}; t_ySBC{i+2}],1);
%         counter= counter +1;
%     end
%     
%       counter= 1;
%     for i =1:4:size(t_yONPar,1)-2
%         tcONPar(counter,:) = mean([t_yONPar{i}; t_yONPar{i+1}; t_yONPar{i+2}],1);
%         counter= counter +1;
%     end
%     
%     
%       counter= 1;
%     for i =1:4:size(t_yOFFPar,1)-2
%         tcOFFPar(counter,:) = mean([t_yOFFPar{i}; t_yOFFPar{i+1}; t_yOFFPar{i+2}],1);
%         counter= counter +1;
%     end
%     
%     
%       counter= 1;
% 
%     for i =1:4:size(t_yONMid,1)-2
%         tcONMid(counter,:) = mean([t_yONMid{i}; t_yONMid{i+1}; t_yONMid{i+2}],1);
%         counter= counter +1;
%     end
%     
%       counter= 1;
%     for i =1:4:size(t_yOFFMid,1)-2
%         tcOFFMid(counter,:) = mean([t_yOFFMid{i}; t_yOFFMid{i+1}; t_yOFFMid{i+2}],1);
%         counter= counter +1;
%     end
%     
%       counter= 1;
%     for i =1:4:size(t_yONLar,1)-2
%         tcONLar(counter,:) = mean([t_yONLar{i}; t_yONLar{i+1}; t_yONLar{i+2}],1);
%         counter= counter +1;
%     end    
% %     tcOFFPar = tcOFFPar*scale_OFF_by;
% %         tcOFFMid = tcOFFMid*scale_OFF_by;
% %         
% %     normalized_time = time_green./repmat(spike_counts,1,30);
% tc = [tcSBC; tcONPar; tcOFFPar; tcONMid; tcOFFMid];
%     meanAcrossTime = mean(tc,1);  
%     normalized_time = tc - repmat(meanAcrossTime, size(tc,1),1);
%        
%        
% % Plot in the their row of the output figure
% %         subplot(3,2+N, 2+N+2+N +2 +i)
% %         
% %         % Plot in ms and only for positive time
% %         plot(1000*time(time>0), ccf(i, time>0))
% %         xlabel('time (ms)')
%     square = normalized_time'*normalized_time;
%     [V,D,W] = eig(square);  
%     [max_eig, ind] = max(diag(D))
%     
%     cell_indices = [cell_indicesSBC,cell_indicesONPar, cell_indicesOFFPar, cell_indicesONMid, cell_indicesOFFMid];
%    
%     TF = zeros(size(tc,1),1);
%     for i = 1:length(cell_indices)
%         TF(i) = normalized_time(i,:)*V(:,ind);
%     end
%     
%     
%     TFSBC = TF(1:size(tcSBC,1));
%     TFONPar = TF(size(tcSBC,1)+1:size(tcSBC,1) + size(tcONPar,1));
%         TFOFFPar = TF(size(tcSBC,1) + size(tcONPar,1) + 1:size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1));
%         TFONMid = TF(size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1)+1:size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1));
%         TFOFFMid = TF(1+ size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1): end);
%         
%       onscale = max(TFONMid);  
%             offscale = max(-TFOFFPar);  
% 
%       TFONMid = TFONMid/onscale;
%       TFONPar = TFONPar/onscale;
%       
%       TFOFFPar = TFOFFPar/offscale;
%       TFOFFMid = TFOFFMid/offscale;
% 
%     area = [areaSBC; areaONPar; areaOFFPar; areaONMid; areaOFFMid];
%     cmap = hsv(6);
%     figure; 
%     plot(areaSBC, TFSBC, '.', 'Color', cmap(1,:),'MarkerSize',25)
%     hold on
%         plot(areaONPar, TFONPar,  '.k', 'MarkerSize', 25)
%         plot(areaONMid, TFONMid,  '.', 'Color', cmap(3,:),'MarkerSize', 25)
%         plot(areaOFFPar, TFOFFPar,  '.', 'Color', cmap(4,:),'MarkerSize', 25)
%         plot(areaOFFMid, TFOFFMid,  '.', 'Color', cmap(5,:),'MarkerSize', 25)
%     set(gca, 'ylim', [-1.2 1.2])
%     set(gca, 'xtick', [0:50:200])
% xlabel('RF diameter (\mum)')
% ylabel('TF_1')
% set(gcf, 'Position', [100,100, 500,750])
%     hgexport(fig, [filepath,'big 5'])
% 
%         
%         %%%% with ON Smooth
% 
%         
%         tc = [tcSBC; tcONPar; tcOFFPar; tcONMid; tcOFFMid; tcONLar];
%     meanAcrossTime = mean(tc,1);  
%     normalized_time = tc - repmat(meanAcrossTime, size(tc,1),1);
%        
%        
% % Plot in the their row of the output figure
% %         subplot(3,2+N, 2+N+2+N +2 +i)
% %         
% %         % Plot in ms and only for positive time
% %         plot(1000*time(time>0), ccf(i, time>0))
% %         xlabel('time (ms)')
%     square = normalized_time'*normalized_time;
%     [V,D,W] = eig(square);  
%     [max_eig, ind] = max(diag(D))
%     
%     cell_indices = [cell_indicesSBC,cell_indicesONPar, cell_indicesOFFPar, cell_indicesONMid, cell_indicesOFFMid, cell_indicesONLar];
%    
%     TF = zeros(size(tc,1),1);
%     for i = 1:length(cell_indices)
%         TF(i) = normalized_time(i,:)*V(:,ind);
%     end
%     
%     
%     TFSBC = TF(1:size(tcSBC,1));
%     TFONPar = TF(size(tcSBC,1)+1:size(tcSBC,1) + size(tcONPar,1));
%         TFOFFPar = TF(size(tcSBC,1) + size(tcONPar,1) + 1:size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1));
%         TFONMid = TF(size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1)+1:size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1));
%         TFOFFMid = TF(1+ size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1):  size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1) + size(tcOFFMid,1));
%                 TFONLar = TF(1+ size(tcSBC,1) + size(tcONPar,1) + size(tcOFFPar,1) + size(tcONMid,1) + size(tcOFFMid,1): end);
% 
%       onscale = max(abs(TFONMid));  
%             offscale = max(abs(TFOFFPar));  
% 
%       TFONMid = TFONMid/abs(onscale);
%       TFONPar = TFONPar/abs(onscale);
%             TFONLar = TFONLar/abs(onscale);
% 
%       TFOFFPar = TFOFFPar/abs(offscale);
%       TFOFFMid = TFOFFMid/abs(offscale);
% 
%     area = [areaSBC; areaONPar; areaOFFPar; areaONMid; areaOFFMid; areaONLar];
%     figure; 
%  plot(areaSBC, TFSBC, '.', 'Color', cmap(1,:),'MarkerSize', 20)
%     hold on
%         plot(areaONPar, TFONPar,  '.k', 'MarkerSize', 20)
%         plot(areaONMid, TFONMid,  '.', 'Color', cmap(3,:),'MarkerSize', 20)
%         plot(areaOFFPar, TFOFFPar,  '.', 'Color', cmap(4,:),'MarkerSize', 20)
%         plot(areaOFFMid, TFOFFMid,  '.', 'Color', cmap(5,:),'MarkerSize', 20)
%         plot(areaONLar, TFONLar,  '.', 'Color', cmap(6,:),'MarkerSize', 20)
%     set(gca, 'xtick', [0:100:550])
% 
% %     set(gca, 'ydir', 'reverse')
%     set(gca, 'ylim', [-1.2 1.2])
%         set(gca, 'xlim', [0 550])
% set(gcf, 'Position', [100,100, 500,750])
% 
% xlabel('RF diameter (\mum)')
% ylabel('TF_1')
% 
% hgexport(fig, [filepath,'big 6'])


% %% ON parasol with ON smooth
% 
% 
%     run_opts.date='2005-04-26-0/'; % one slash at the end
% run_opts.concatname='data009-dum'; % Name (or modified name) of run, no slashes
% run_opts.concatname_parasol='data009-dum'; % Name (or modified name) of run, no slashes
% 
% % Sometimes the data has two versions of the concate name
% % run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
% run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data009/data009'];
% 
% % Full path to movie xml
% run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-1-0.48-33333.xml';
% 
% 
% run_opts.save_location_root = '/Users/colleen/Desktop/Berstein Talk/Mosaics/';
% % Number of frames to use for generator signal as well as number of frames
% % of the timecourse to display
% run_opts.num_frames = 25;
% 
% % Number of bins to use for the nonlinearity graph
% run_opts.num_bins = 10;
% 
% % How much padding to use for zooming in on STAs
% params.padding = 7;
% 
% % Cell specification can be one cell type or multiple in a cell array.
% % Use the same spelling/capitalization as the vision params file
% % OFF large-2 in vision = OFF large 2
% cell_specification = {'OFF parasol', 'OFF large 1' ,'OFF large 2' };
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where to save
% run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];
% 
% % Used for labeling plot
% run_opts.cell_type = cell_specification;
% 
% % Load Data
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
% opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
% datarun=load_data(datarun,opt);
% 
% filepath=[sprintf('%s%s%s.pdf',[run_opts.save_location_root,'/',run_opts.date,'/',run_opts.concatname, '/'])];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end
% 
%     
%     
%     
%     %% Plot mosaics on top of each other
%     
%     
%     
%    run_opts.load_location = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname];
% run_opts.load_location_parasol = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname_parasol];
% 
% % load fitted datadata014
% s=dir(fullfile(run_opts.load_location,'*.mat'));
% 
%     for cell_spec = 1:size(cell_specification,2)
%         cells = cell_specification{cell_spec};
%         count = 1;
%         for i = 1:size(s,1)
%             period = strfind(s(i).name, '.');
%             if sum(strcmp(s(i).name(1:period-1), cells))>0
%                 fitted_cells{cell_spec}(count) =  load([run_opts.load_location, '/', s(i).name]);
%                 count = count+1;
%             end
%             
%             
%         end
%     end
%     
% 
%     
%     
%     
%       %%%%%%%%%%%%%%%%%%%%%%% More than one Cell Type %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Generate an output figure that will have the overlaid mosaics and
%     % timecourses of the inputted cell types
%     fig = figure;
%     set(fig, 'Position', [0 0 700 700])
%     set(0,'DefaultAxesFontSize',10)
%     
%     % Will have the text for the legend
%     Legend= cell(size(cell_specification,2),1);
%     
%     for i = 1:size(cell_specification,2)
%         
%         % For each cell type
%         cell_spec = cell_specification{i};
%         
%         % color to plot this cell type in
%         cmap = distinguishable_colors(size(cell_specification,2));
%         cmap(1,:) = [0,0,0];
%         
%         % Now many cells are in this type
%         [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_spec);
%         cell_ids=datarun.cell_ids(cell_indices);
%         N = length(cell_ids);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         % Add the mosaic to the output figure
%         figure(fig) % make the output figure the current gcf
%         hold on;
%         subplot(3,2, [1,2,3,4])
%         % Get the axes handle for a figure where all the rf fits are plotted
%         rf_c = get(fitted_cells{i}.prop.mosaic_axes, 'Children');
%         rf_x = get(rf_c,'Xdata');
%         rf_y = get(rf_c,'Ydata');
%         
%         
%         
%         
%         
%         
%         % type of rf_x and rf_y is only a cell array when there is more than
%         % one cell in the type
%         if N>1
%             % Plot each cell type in a different color
%             if i ~= 1
%                 h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:),  'linewidth', 2);
%                 
%                 
%             else
%                 h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:));
%                 
%             end
%         else
%             if i ~= 1
%                 h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:), 'linewidth', 2);
%                 
%                 
%             else
%                 h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:));
%                 
%             end
%         end
%         
%         % Add the name of the cell type to the legend cell array
%         Legend{i} = cell_spec;
%         
%         % When in the last iteration of the loop
%         if i == size(cell_specification,2)
%             for j = 1:size(cell_specification,2)
%                 % Get the line property of the first cell of each cell type
%                 % so that the colors in the legend are correct
%                 inter(j) = h{j}(1);
%             end
%             
%             % Specify the data and the string for each cell type
%             legend_handle = legend([inter] , Legend);
%         end
%         hold on
%         axis equal
%         % Set the axis of the mosaic to what the stimulus size was
%         %     axis([bounds])
% %         axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
%         set(gca, 'ydir', 'reverse')
%         
%         % make the output figure invisible, close the figure that came from
%         % plot_rf_summaries, and then turn the output figure back on
%         % Pretends tons of figure from being open at the end of the script
%         set(fig, 'HandleVisibility', 'off');
%         close all;
%         set(fig, 'HandleVisibility', 'on');
%         
%     end
%     
%     
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Timecourse of cell type i %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     for i = 1:size(cell_specification,2)
%         % For each cell type, get the overlaid timecourse
%         cell_spec = cell_specification{i};
%         params.normalize = 1;
%         %         t = plot_all_timecourses(datarun, cell_spec, params, run_opts);
%         % Make the output figure the gcf
%         figure(fig)
%         subplot(3,2, [5,6]) % Plot 1x2 subplot
%         
%         t_c = get(fitted_cells{i}.prop.t_axes, 'Children');
%         
%         %         t_c = get(t, 'Children');
%         %         xtick = get(t, 'xtick');
%         %         xticklabel = get(t, 'xticklabel');
%         
%         t_x = get(t_c,'Xdata');
%         t_y = get(t_c,'Ydata');
%         time = round(fliplr(0:-datarun.stimulus.interval*8.33:-(run_opts.num_frames-1)*datarun.stimulus.interval*8.33));
%         for j = 1:size(t_x,1)
%             
%             if get(t_c(j), 'color') == [0,1,0]
%                 try
%                     H{i} = plot(time, t_y{j}(end-run_opts.num_frames+1:end), 'color', [cmap(i,:), 0.5]);
%                     hold on
%                 catch
%                     H{i} =plot(time, t_y(end-run_opts.num_frames+1:end), ['color',cmap(i,:), 0.5]);
%                     hold on
%                 end
%                 
%             end
%             
%         end
%         
%         %         H{i} =  plot(time, fitted_cells.prop.avg_timecourse(end-run_opts.num_frames+1:end), 'k')
%         set(gca, 'xtick', time(1:2:end))
%         set(gca, 'xticklabel', time(1:2:end));
%         
%         axis tight
%         
%         hold on
%         
%         
%         Legend{i} = cell_spec;
%         
%         % if you are on the last iteration of the loop
%         if i == size(cell_specification,2)
%             % Select the line axes properites of the first cell of each
%             % type
%             for j = 1:size(cell_specification,2)
%                 inter(j) = H{j}(1);
%             end
%             
%             % Set the line properties and string properites of the legend
%             legend_handle = legend([inter] , Legend, 'Location', 'Southwest');
%             
%         end
%     end
%     plot(get(gca, 'xlim'), zeros(1,2), 'k-')
%     xlabel('time (ms)')
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Saving the output figure %%%%%%%%%%%%%%%%%%%%%%
%     % Name the figure the overall cell type category (ie ON large when the
%     % cell_specification is {'ON large 1' 'ON large 2'}
%     spaces = strfind(cell_specification{2}, ' ');
%     title_phrase = cell_specification{2}(1:(spaces(2)-1));
%     suptitle({[run_opts.date, run_opts.concatname]; title_phrase})
%     
%     % Make the directory to save in if it doesn't exist
%     if ~exist(run_opts.filepath)
%         mkdir(run_opts.filepath)
%     end
%     
%     % Save the figure to pdf
%     export_fig([run_opts.filepath, [title_phrase, ' Summary']], '-pdf')
%     
%     
%     %%  per cell type
%     cell_specification = {'OFF parasol'};
% 
%     
%  
%     % To figure out how many cells are in the cell type
%     [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
%     cell_ids=datarun.cell_ids(cell_indices);
%     N = length(cell_ids);
%     
%     
% %        Final figure where each component will be plotted
%     fig = figure;
%     % Size depends on number of cell you are plotting
%     set(fig, 'Position', [0 0 1000 300])
%     set(0,'DefaultAxesFontSize',10)
%     
%   % Get the axes handle for a figure where all the rf fits are plotted
%     rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
%     
%         % Add the mosaic to the output figure
%     figure(fig) % make the output figure the current gcf
%     if N > 5
%         N = 5;
%     end
%     
%     subplot(3,2+N, [1,2,2+N+1,2+N+2]) % Top left 2x2 subplot
%     
%     % get the data plotted in plot_rf_summaries
%     rf_c = get(rf, 'Children');
%     rf_x = get(rf_c,'Xdata');
%     rf_y = get(rf_c,'Ydata');
%     
%     % type of rf_x and rf_y is only a cell array when there is more than
%     % one cell in the type
%     if N>1
%         plot(cell2mat(rf_x)', cell2mat(rf_y)', 'k')
%     else
%         plot(rf_x, rf_y, 'k')
%     end
%     
%     
%     %     if strcmp(datarun.stimulus.independent, 'nil')
%     %         colormap gray
%     %     end
%     
%     axis equal
%     % Set the axis of the mosaic to what the stimulus size was
% %     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
% %     axis([ 20 50 0 25])
% 
%     set(gca, 'ydir', 'reverse')
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_rf_summaries, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
% %     set(fig, 'HandleVisibility', 'off');
% %     close all;
% %     set(fig, 'HandleVisibility', 'on');
% 
% 
%  %%%%%%%%%%%%%%%%%%%%%%% Overlaid Timecourses %%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Turns the axes handle to a plot of the overlaid timecourses
%     params.normalize = 1;
%     t = plot_all_timecourses(datarun, cell_specification, params, run_opts);
%     
%     % Make the output figure the gcf
%     figure(fig)
%     subplot(3,2+N, [2*(2+N)+1,2*(2+N)+2]) % Lower left 2x2 plot
%     
%     % Extract the lines from the plot from plot_all_timecourses
%     t_c = get(t, 'Children');
%     xtick = get(t, 'xtick');
%     xticklabel = get(t, 'xticklabel');
%     
%     t_x = get(t_c,'Xdata'); % obtain the XData
%     t_y = get(t_c,'Ydata'); % obtain the YData
%     
%     % Plot all the blue channels together
%     index= zeros(size(t_x, 1),1);
%     index(1:4:end) = 1; % blue information is 1,4,7 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1),:), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [0 0 1 0.55])
%     
%     hold on
%     % Plot all the red channels together
%     
%     index(1:4:end) = 0;
%     index(3:4:end) = 1; % red information is 3,6,9 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [1 0 0 0.55])
%     
%     % Plot all the green channels together
%     % Plot green last to it shows in front of the other channels
%     index(3:4:end) = 0;
%     index(2:4:end) = 1;% green information is 2,5,8 etc
%     x = cell2mat(t_x(logical(index)))';
%     y = cell2mat(t_y(logical(index)))';
%     % Only plot a certain num_frames
%     
%     plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:),'Color', [0 1 0 0.55])
%     
%     % Plot a line at 0 to help see the zero crossing
%     plot(get(gca, 'xlim'), zeros(1,2), 'k-')
%     
%     axis tight
%     xlabel('time (ms)')
%     
%     % reduce axis labels to just num_frames
%     set(gca, 'xtick', xtick((size(x,1)-run_opts.num_frames):2:size(x,1)));
%     set(gca, 'xticklabel', xticklabel((size(x,1)-run_opts.num_frames):2:size(x,1),:));
%     
%     % make the output figure invisible, close the figure that came from
%     % plot_all_timecourses, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     
%     
%     %%%%%% Individual STAs%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  % Get cell array of axes handles to the plots of the stas zoomed in
%     % BUG: if the cell is at the edge of the array the location of the
%     % scale bar is wrong and the size looks funky
%     p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
%     if N > 5 
%         N = 5;
%     end
%     
%     % Add individual stas to the output figure
%     figure(fig)
%     
%     for i = 1:length(cell_ids)
%         subplot(3,2+N,2+i) % First row subplots 3...12
%         
%         % Get information about the ith cell's STA
%         sta_c = get(p_sta{i}, 'Children');
%         sta_image = get(sta_c(3),'Cdata');
%         % Plot the image before the fit so it is in the background
%         imagesc(sta_image);
%         hold on;
%         
%         % Plot the RF fit
%         sta_1x = get(sta_c(1),'Xdata');
%         sta_1y = get(sta_c(1),'Ydata');
%         plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)
%         
%         % Plot the scale bar (set to 15 pixels)
%         sta_2x = get(sta_c(2),'Xdata');
%         sta_2y = get(sta_c(2),'Ydata');
%         hold on;
%         plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
%         axis([get(p_sta{i}, 'xlim') get(p_sta{i}, 'ylim')])
%         axis square
%         axis off
%         
%         % Don't plot more than 10 cells (not enough space)
%         if i == 5
%             break;
%         end
%         
%     end
%     
%         % make the output figure invisible, close the figure that came from
%     % plot_zoomed_STAs, and then turn the output figure back on
%     % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');
%     





%% ON parasol with multiple ON large


    run_opts.date='2007-05-01-2/'; % one slash at the end
run_opts.concatname='data000'; % Name (or modified name) of run, no slashes
run_opts.concatname_parasol='data000'; % Name (or modified name) of run, no slashes

% Sometimes the data has two versions of the concate name
% run_opts.file_name = [run_opts.date, '/', 'data006-nwpca', '/',  'data006/data006'];
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  'data000'];

% Full path to movie xml
run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-4-0.48-11111.xml';


run_opts.save_location_root = '/Users/colleen/Desktop/Berstein Talk/Mosaics/';
% Number of frames to use for generator signal as well as number of frames
% of the timecourse to display
run_opts.num_frames = 25;

% Number of bins to use for the nonlinearity graph
run_opts.num_bins = 10;

% How much padding to use for zooming in on STAs
params.padding = 7;

% Cell specification can be one cell type or multiple in a cell array.
% Use the same spelling/capitalization as the vision params file
% OFF large-2 in vision = OFF large 2
cell_specification = {'ON parasol', 'ON large 1' ,'ON large 2' };



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

filepath=[sprintf('%s%s%s.pdf',[run_opts.save_location_root,'/',run_opts.date,'/',run_opts.concatname, '/'])];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

    
    
    
    %% Plot mosaics on top of each other
    
    
    
   run_opts.load_location = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname];
run_opts.load_location_parasol = ['/Users/colleen/Desktop/Cell Characteristics/', run_opts.date, '/', run_opts.concatname_parasol];

% load fitted datadata014
s=dir(fullfile(run_opts.load_location,'*.mat'));

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
        cmap(1,:) = [0.1,0.1,0.1];
        
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
            if i ~= 1
                h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:),  'linewidth', 2);
                
                
            else
                h{i} =  plot(cell2mat(rf_x)', cell2mat(rf_y)', 'Color', cmap(i,:));
                
            end
        else
            if i ~= 1
                h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:), 'linewidth', 2);
                
                
            else
                h{i} =  plot(rf_x, rf_y, 'Color', cmap(i,:));
                
            end
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
%         axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
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
                    H{i} = plot(time, t_y{j}(end-run_opts.num_frames+1:end), 'color', [cmap(i,:), 0.5]);
                    hold on
                catch
                    H{i} =plot(time, t_y(end-run_opts.num_frames+1:end), ['color',cmap(i,:), 0.5]);
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
    
    
    %%  per cell type
    cell_specification = {'ON large 2'};

    
 
    % To figure out how many cells are in the cell type
    [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
    cell_ids=datarun.cell_ids(cell_indices);
    N = length(cell_ids);
    
    
%        Final figure where each component will be plotted
    fig = figure;
    % Size depends on number of cell you are plotting
    set(fig, 'Position', [0 0 1000 300])
    set(0,'DefaultAxesFontSize',10)
    
  % Get the axes handle for a figure where all the rf fits are plotted
    rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
    
        % Add the mosaic to the output figure
    figure(fig) % make the output figure the current gcf
    if N > 5
        N = 5;
    end
    
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
            colormap(gray)

    
    %     if strcmp(datarun.stimulus.independent, 'nil')
    %         colormap gray
    %     end
    
    axis equal
    % Set the axis of the mosaic to what the stimulus size was
%     axis([ 0 datarun.stimulus.field_width 0 datarun.stimulus.field_height])
%     axis([ 20 50 0 25])

    set(gca, 'ydir', 'reverse')
    
    % make the output figure invisible, close the figure that came from
    % plot_rf_summaries, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
%     set(fig, 'HandleVisibility', 'off');
%     close all;
%     set(fig, 'HandleVisibility', 'on');


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
    plot(x((size(x,1)-run_opts.num_frames):size(x,1),:), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [0 0 0 0.55])
    
    hold on
    % Plot all the red channels together
    
    index(1:4:end) = 0;
    index(3:4:end) = 1; % red information is 3,6,9 etc
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    % Only plot a certain num_frames
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:), 'Color', [0 0 0 0.55])
    
    % Plot all the green channels together
    % Plot green last to it shows in front of the other channels
    index(3:4:end) = 0;
    index(2:4:end) = 1;% green information is 2,5,8 etc
    x = cell2mat(t_x(logical(index)))';
    y = cell2mat(t_y(logical(index)))';
    % Only plot a certain num_frames
    
    plot(x((size(x,1)-run_opts.num_frames):size(x,1), :), y((size(x,1)-run_opts.num_frames):size(x,1),:),'Color', [0 0 0 0.55])
    
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
    
    
    %%%%%% Individual STAs%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Get cell array of axes handles to the plots of the stas zoomed in
    % BUG: if the cell is at the edge of the array the location of the
    % scale bar is wrong and the size looks funky
    p_sta = plot_zoomed_STAs(datarun, cell_specification, params);
    if N > 5 
        N = 5;
    end
    
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
        if i == 6
            break;
        end
        
    end
    
        % make the output figure invisible, close the figure that came from
    % plot_zoomed_STAs, and then turn the output figure back on
    % Pretends tons of figure from being open at the end of the script
    set(fig, 'HandleVisibility', 'off');
    close all;
    set(fig, 'HandleVisibility', 'on');
        export_fig([run_opts.filepath, cell_specification{1}], '-pdf')

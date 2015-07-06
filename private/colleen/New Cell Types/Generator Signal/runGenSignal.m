% runGenSignal

clear


run_opts.date='2010-09-24-0/';
run_opts.concatname='data001-nwpca-cr';
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-80x60.xml';

run_opts.cell_type = {'OFF large 2'};
run_opts.num_frames = 19;
run_opts.num_bins = 10; % Number of bins to use for the generator vs spikes graph

% cell_specification_large = [21]; % visionID
cell_specification = {'OFF large 2'}; % visionID

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_opts.filepath= ['/Users/colleen/Desktop/Generator Signal/', run_opts.date, '/', run_opts.concatname, '/'];


% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames =1:30;% have to input as a vector list of frames, not the number of frames total
datarun=load_data(datarun,opt);
% datarun = load_params(datarun,'verbose',1);
% datarun = load_neurons(datarun);
% datarun = load_sta(datarun);
% datarun = set_polarities(datarun);


% Mosaic
rf = plot_rf_summaries(datarun, cell_specification, 'foa', 0);
% plot_time_course(datarun, [21 33 76], 'figure', 0)

% Timecourses
params.normalize = 1;
t = plot_all_timecourses(datarun, cell_specification, params);

% STAs 
params.padding = 10;
p_sta = plot_zoomed_STAs(datarun, cell_specification, params)

% Gen Signal
gs = computeAndPlot_genSignal(datarun, cell_specification, run_opts );

%% ACF
[cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
 cell_ids=datarun.cell_ids(cell_indices);
 for i = 1:length(cell_ids)
     cell_id = cell_ids(i);
    [time, ccf(i,:)] = plot_ccf(datarun, cell_id);
 end
 avg_ccf = mean(ccf,1);
 figure
plot(1000*time(time>0), avg_ccf(time>0))

%% Mosaic
rf_c = get(rf, 'Children');
rf_x = get(rf_c,'Xdata'); % obtain the XData 
rf_y = get(rf_c,'Ydata'); % obtain the YData 
 figure; plot(cell2mat(rf_x)', cell2mat(rf_y)')
 axis([ 0 80 0 60])
 set(gca, 'ydir', 'reverse')
 
 
 %% Timecourse
 t_c = get(t, 'Children');
t_x = get(t_c,'Xdata'); % obtain the XData 
t_y = get(t_c,'Ydata'); % obtain the YData 

index= ones(size(t_x, 1),1);
index(4:4:end) = 0;
 figure; plot(cell2mat(t_x(logical(index)))', cell2mat(t_y(logical(index)))')
 
 
 %% STAs
 
 sta_c = get(p_sta, 'Children');
   figure;
   hold on;
  sta_image = get(sta_c(3),'Cdata');
imagesc(sta_image);

 sta_1x = get(sta_c(1),'Xdata');
 sta_1y = get(sta_c(1),'Ydata');
 plot(sta_1x, sta_1y,  '-k', 'LineWidth', 2)

   
 sta_2x = get(sta_c(2),'Xdata');
 sta_2y = get(sta_c(2),'Ydata');
hold on;
plot(sta_2x, sta_2y, 'Color','k', 'LineWidth',1)
axis([get(p_sta, 'xlim') get(p_sta, 'ylim')])
 
%% Generator Signal

gs_c = get(gs, 'Children');
gs_x = get(gs_c,'Xdata'); % obtain the XData 
gs_y = get(gs_c,'Ydata'); % obtain the YData 
 figure; plot(gs_x, gs_y)

 

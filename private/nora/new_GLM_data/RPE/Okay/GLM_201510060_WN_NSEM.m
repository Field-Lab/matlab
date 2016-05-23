% RPE, WN and NSEM
   %{
  clear
Analysis_Path = '/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/';
class = 'data000';
NSEM_runs = {'data002', 'data007', 'data008', 'data009', 'data012', 'data015'};
NSEM_movie = {'A', 'B', 'D', 'E', 'F', 'C'};
WN_runs = {'data003', 'data010', 'data013'};
WN_seeds = [11111 22222 33333];
WN_fit_frames = [3600 6000 3600];
extra_string = '-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015';
datarun_class = load_data([Analysis_Path class extra_string '/' class extra_string], struct('load_neurons', 0, 'load_params', 1));
mkdir('/Volumes/Lab/Users/Nora/GLMFits/RPE/201510060/Midget/');
mkdir('/Volumes/Lab/Users/Nora/GLMFits/RPE/201510060/Midget/WN/');
%mkdir('/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/');
 %}
%% WN
%%{
reps = 60;
idx = 1:reps;
cell_type = {'Off Midget Clean'};
cells = get_cell_ids(datarun_class, cell_type{1}); % cell ids to fit
% on midget
%cells = cells([13 47 51 53 56 82 90 92 95]);
% datarun_class = load_sta(datarun_class, 'load_sta', cells);
n_cells = length(cells);

%%
  %{
run = 1;
datarun = load_data([Analysis_Path WN_runs{run} extra_string '/' WN_runs{run} extra_string], struct('load_neurons', 1, 'load_params', 1));
prepped_data = interleaved_data_prep(datarun, [WN_fit_frames(run) 1200], reps,'cell_spec', cells,'visual_check', 0);
%}
for run = 1:3
    datarun = load_data([Analysis_Path WN_runs{run} extra_string '/' WN_runs{run} extra_string], struct('load_neurons', 1, 'load_params', 1));
 prepped_data = interleaved_data_prep(datarun, [WN_fit_frames(run) 1200], reps,'cell_spec', cells,'visual_check', 0);
    fitspikes_cell(idx, :) = prepped_data.fitspikes;
    % fitmovie_cell(idx) = prepped_data.fitmovie;
    testspikes_cell(idx, :) = prepped_data.testspikes;
    % testmovie = prepped_data.testmovie;
    monitor_refresh(run) = 100/median(diff(datarun.triggers));
    idx = idx+reps;
end

%%
monitor_refresh = mean(monitor_refresh);

 for i=1:n_cells
fitspikes{i} = [];
idx = (1:60);
block_lengths = [3600 6000 3600]/monitor_refresh;
offsets = [0 3600*60 3600*60+6000*60]/monitor_refresh;
for run = 1:3
    temp = concat_spikes(fitspikes_cell(idx,i), block_lengths(run));
    temp = temp+offsets(run);
    fitspikes{i} = [fitspikes{i}; temp];
    idx = idx+60;
end
end
fitmovie = concat_movie(fitmovie_cell);

%%
for cell = 1:n_cells
cell_idx = get_cell_indices(datarun_class, cells(cell));
vision_center = datarun_class.vision.sta_fits{cell_idx}.mean;
center = round([40-vision_center(2), vision_center(1)]);
    fittedGLM = glm_fit(fitspikes{cell}, fitmovie, center, 'monitor_refresh', monitor_refresh);
    fittedGLM.xvalperformance = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes_cell(:,cell));
    temp = corrcoef(conv(sum(fittedGLM.xvalperformance.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xvalperformance.rasters.recorded), gausswin(100)));
    fittedGLM.xvalperformance.corr = temp(2,1);
    close all
    plotfilters(fittedGLM)
    exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/RPE/201510060/Midget/WN/OffMid_' num2str(cells(cell)) '_filters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    plotrasters(fittedGLM.xvalperformance, fittedGLM)
    exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/RPE/201510060/Midget/WN/OffMid_' num2str(cells(cell)) '_rasters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    save(['/Volumes/Lab/Users/Nora/GLMFits/RPE/201510060/Midget/WN/OffMid_' num2str(cells(cell)) '.mat'], 'fittedGLM');
end
%}

%% NSEM
%{
clear fitspikes_cell
clear testspikes_cell
reps = 50; 
idx = 1:reps;
cell_type = {'Off Midget'};
cells = get_cell_ids(datarun_class, cell_type{1}); % cell ids to fit
datarun_class = load_sta(datarun_class, 'load_sta', cells);
n_cells = length(cells);

for run = 1
    datarun = load_data([Analysis_Path NSEM_runs{run} extra_string '/' NSEM_runs{run} extra_string], struct('load_neurons', 1, 'load_params', 1));
    prepped_data = interleaved_data_prep(datarun, [7200 3600], reps,'cell_spec', cells,'visual_check', 1);
    fitspikes_cell(idx, :) = prepped_data.fitspikes;
    monitor_refresh(run) = 100/median(diff(datarun.triggers));
    load(['/Volumes/Lab/Users/Nora/Stimulus/NSEM_Movies/downsampledNSbrownian_' NSEM_movie{run} '_3000.mat']);
    if run == 1
        testspikes_cell(idx, :) = prepped_data.testspikes;
        testmovie = fitmovie(:,:,1:3600);
    end
    fitmovie = fitmovie(:,:,3601:end);
    idx = idx+reps;
end
monitor_refresh = mean(monitor_refresh);

%%
for i=1:n_cells
fitspikes{i} = [];
idx = (1:reps);
block_lengths = 7200/monitor_refresh;
offsets = [0 3600*60 3600*60+6000*60]/monitor_refresh;
for run = 1
    temp = concat_spikes(fitspikes_cell(idx,i), block_lengths);
    temp = temp+offsets(run);
    fitspikes{i} = [fitspikes{i}; temp];
    idx = idx+60;
end
end
close all
for cell = 1:n_cells
%[STA, center] = STA_Test(fitspikes{cell}, fitmovie, 1, 1/monitor_refresh);
cell_idx = get_cell_indices(datarun_class, cells(cell));
vision_center = datarun_class.vision.sta_fits{cell_idx}.mean;
center = round([40-vision_center(2), vision_center(1)]);
fittedGLM = glm_fit(fitspikes{cell}, fitmovie, center, 'monitor_refresh', monitor_refresh);
fittedGLM.xval = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes_cell(:,cell));
temp = corrcoef(conv(sum(fittedGLM.xval.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xval.rasters.recorded), gausswin(100)));
fittedGLM.xval.corr = temp(2,1);
close all
plotfilters(fittedGLM)
exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/OffMid_' num2str(cells(cell)) '_filters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
close all
plotrasters(fittedGLM.xval, fittedGLM)
exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/OffMid_' num2str(cells(cell)) '_rasters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
close all
save(['/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/OffMid_' num2str(cells(cell)) '.mat'], 'fittedGLM');

end
%}

%%
%{
files = dir('/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/OffMid_*.mat');
n_files = length(files);
for i=1:n_files
    cell_name = files(i).name(7:10);
    if strcmp(cell_name(end), '.')
        cell_name = cell_name(1:(end-1));
    end
    load(['/Volumes/Lab/Users/Nora/GLMFits/201510060/NSEM/' files(i).name]);
    temp = corrcoef(conv(sum(fittedGLM.xval.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xval.rasters.recorded), gausswin(100)));
    MSE_NSEM(i) = temp(2,1);
end
%}

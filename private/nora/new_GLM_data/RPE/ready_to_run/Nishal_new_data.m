clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
cells = get_cell_ids(datarun_class, 'Off Parasol'); % cell ids to fit

%% WN 
%%{
test_data = 'data027';
fit_data = 'data026';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 22222);
FG = glm_fit_from_WN(cells(1), [Analysis_Path '/' fit_data '/' fit_data], 'BW-8-1-0.48-11111', 'testmovie', repeats.testmovie, 'testspikes', repeats.testspikes, 'monitor_refresh', monitor_refresh, 'stim_length', 1800);
disp('done with WN glm fit')
%}
%% NSEM
%%{
test_data = 'data022';
fit_data = 'data025';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30,'cell_spec', cells,'visual_check', 0);
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
monitor_refresh = 100/median(diff(fit_datarun.triggers));
load('/Volumes/Lab/Users/Nora/Stimulus/NSEM_Movies/downsampledNSinterval.mat');
testmovie = fitmovie(:,:,1:1100);
fitmovie = fitmovie(:,:,1101:end);

for i = 1:length(cells)
    disp(i)
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
    fitspikes = fitspikes(fitspikes>1100/monitor_refresh)-1100/monitor_refresh;
    close all 
    [STA, ~] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
    center(2) = datarun_class.vision.sta_fits{master_idx}.mean(1);
    center(1) = 40 - datarun_class.vision.sta_fits{master_idx}.mean(2);
    center = round(center);
    close all
    fittedGLM     = glm_fit(fitspikes, fitmovie, center, 'monitor_refresh', monitor_refresh);
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    plotfilters(fittedGLM);
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
end

%}

clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_test/2016-02-17-1';
monitor_refresh = 119.5;

%% only change things here
cells = get_cell_ids(datarun_class, 'On Parasol'); % cell ids to fit

%% NSEM
test_data = 'data022';
fit_data = 'data025';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30,'cell_spec', cells,'visual_check', 0);
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat')
%stim_length = (3600-10)*convergence;
%fitmovie = fitmovie(:,:,1:ceil(stim_length*monitor_refresh));
testmovie = fitmovie(:,:,1:1200);

%%

for i = 1%:3%length(cells)
    
    disp(i)
    glm_cellinfo.cid           = cells(i);
    glm_cellinfo.cell_savename = num2str(cells(i));
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
    %fitspikes = fitspikes(fitspikes < stim_length);
    [STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
    
    
    fittedGLM     = glm_fit(fitspikes, fitmovie, center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
    %{
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    %}
end

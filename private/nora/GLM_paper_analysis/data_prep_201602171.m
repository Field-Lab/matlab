clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits/2016-02-17-1';
monitor_refresh = 119.5;

%% only change things here
cells = get_cell_ids(datarun_class, 'Off Parasol'); % cell ids to fit
convergence = 0.5; % fraction of data to use

%%
if convergence < 1; dsave = [dsave '_Conv_' num2str(convergence)]; end
mkdir(dsave)

%% BW

test_data = 'data027';
fit_data = 'data026';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 22222);
glm_fit_from_WN(cells, [Analysis_Path '/' fit_data '/' fit_data], 'BW-8-1-0.48-11111', 'testmovie', repeats.testmovie, 'testspikes', repeats.testspikes, 'd_save', dsave, 'monitor_refresh', monitor_refresh, 'stim_length', 1800*convergence);
disp('done with WN glm fit')

%% NSEM
test_data = 'data022';
fit_data = 'data025';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30,'cell_spec', cells,'visual_check', 0);
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat')
testmovie = fitmovie(:,:,1:1200);
if convergence < 1
    stim_length = (3600-10)*convergence;
    fitmovie = fitmovie(:,:,1:ceil(stim_length*monitor_refresh));
end


for i = 1:length(cells)
    disp(i)
    glm_cellinfo.cid           = cells(i);
    glm_cellinfo.cell_savename = num2str(cells(i));
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
    fitspikes = fitspikes(fitspikes < stim_length);
    
    eval(sprintf('load %s/%s.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    BPS(i,1) = fittedGLM.xvalperformance.glm_normedbits;
    [STA, center] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
   
    fittedGLM     = glm_fit(fitspikes, fitmovie, fittedGLM.center, 'WN_STA', fittedGLM.STA, 'monitor_refresh', monitor_refresh);
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    eval(sprintf('save %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfilters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
    
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMrasters'], 'Bounds', 'loose', 'Color', 'rgb');
    
    close all

end

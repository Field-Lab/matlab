%%{
clear
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2015-12-18-2/data002-data015';
datarun_class = load_data([Analysis_Path '/data002/data002'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2015-12-18-2/fullscreen/';
mkdir(dsave);
monitor_refresh = 119.5;
%cells = [1925 3181 5794 7246];

% group
cells = [167 1082 2359 3271 5478 6648];
%}

fit_data = 'data003';
test_data = 'data005';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells);
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
load('/Volumes/Data/2015-12-18-2/Visual/Photons_params_timestamps/data003/data003_time_stamps.mat');
TS_trigs = time_stamps{1}(1:100:end)-time_stamps{1}(1)+ fit_datarun.triggers(1);
stim_length = (3100)*convergence;
%fitmov_orig = fitmovie(:,:,1:ceil(stim_length*monitor_refresh));
%fitmovie = uint8(double(fitmovie).*repmat(mask, [1 1 size(fitmovie,3)]) + 64*ones(size(fitmovie)).*(1-repmat(mask, [1 1 size(fitmovie,3)])));
%load('/Volumes/Data/2016-01-05-0/Visual/timestamps/2016-01-05-0/data009_time_stamps.mat');
testmovie = fitmovie(:,:,1:1200);
fitmovie = fitmovie(:,:,1:(3100*120));
%%
for i = 1:6
    disp(i)
    glm_cellinfo.cid           = cells(i);
    glm_cellinfo.cell_savename = num2str(cells(i));
    load(['/Volumes/Lab/Users/Nora/GLMFits_masking/2015-12-18-2/' glm_cellinfo.cell_savename 'NSEM.mat'])
    center = [fittedGLM.center_coord.x_coord fittedGLM.center_coord.y_coord]; clear fittedGLM
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, TS_trigs, 100, monitor_refresh);    
    fitspikes = fitspikes(fitspikes < stim_length);
    %[STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
   
    fittedGLM     = glm_fit(fitspikes, fitmovie,center, 'monitor_refresh', monitor_refresh);%, 'WN_STA', STA);
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    eval(sprintf('save %s/%sNSEM_full.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfilters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
    
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMrasters'], 'Bounds', 'loose', 'Color', 'rgb');
    
    close all

end
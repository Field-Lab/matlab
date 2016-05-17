%%{
clear
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2016-01-05-0/data001-data016';
datarun_class = load_data([Analysis_Path '/data001/data001'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2016-01-05-0';
monitor_refresh = 119.5;
%cells = [1925 3181 5794 7246];

% group
cells = [214 2341 3483 5641];
load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma2/mask_and_filter.mat');
mask = imresize(final_mask,1/4, 'box');
%}

centers = [29 34; 28 57; 11 62; 13 30];

fit_data = 'data009';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
stim_length = (3600)*convergence;
%fitmov_orig = fitmovie(:,:,1:ceil(stim_length*monitor_refresh));
fitmovie = uint8(double(fitmovie).*repmat(mask, [1 1 size(fitmovie,3)]) + 64*ones(size(fitmovie)).*(1-repmat(mask, [1 1 size(fitmovie,3)])));
%load('/Volumes/Data/2016-01-05-0/Visual/timestamps/2016-01-05-0/data009_time_stamps.mat');

%%
for i = 1:4
    disp(i)
    glm_cellinfo.cid           = cells(i);
    glm_cellinfo.cell_savename = num2str(cells(i));
    master_idx         = find(fit_datarun.cell_ids == cells(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun, 100, monitor_refresh);    
    fitspikes = fitspikes(fitspikes < stim_length);
    % eval(sprintf('load %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    % BPS(i,1) = fittedGLM.xvalperformance.glm_normedbits;
    [STA, center] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
   
    fittedGLM     = glm_fit(fitspikes, fitmovie, centers(i,:), 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
    %fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    eval(sprintf('save %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfilters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
    
    %plotrasters(fittedGLM.xvalperformance, fittedGLM);
    %exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMrasters'], 'Bounds', 'loose', 'Color', 'rgb');
    
    close all

end

%%{
clear 
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2016-04-21-1/';
datarun_class = load_data([Analysis_Path 'streamed/data015/data015'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1';
monitor_refresh = 119.5;
cells_fit = [814, 3291, 6050, 4503];
cells_test = [936 3169 6053 4505];
%load('/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma2.mat');
%mask = imresize(mask,1/4, 'box');
%}

fit_data = 'data016';
test_data = 'data018';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
testmovie = fitmovie(:,:,1:1200);
load('/Volumes/Lab/Users/Nora/downsampledNSbrownian_A_1600.mat');
stim_length = (1500)*convergence;
fitmovie = fitmovie(:,:,1:(1500*120));
%fitmovie = uint8(double(fitmovie).*repmat(mask, [1 1 size(fitmovie,3)]) + 64*ones(size(fitmovie)).*(1-repmat(mask, [1 1 size(fitmovie,3)])));
repeats = interleaved_data_prep(test_datarun, 1100, 29, 'cell_spec', cells_test, 'visual_check', 1);

%%
for i = 4
    disp(i)
    glm_cellinfo.cid           = cells_fit(i);
    glm_cellinfo.cell_savename = num2str(cells_fit(i));
    master_idx         = find(fit_datarun.cell_ids == cells_fit(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);    
    fitspikes = fitspikes(fitspikes < stim_length);
    figure;
    [STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
   
    fittedGLM     = glm_fit(fitspikes, fitmovie,center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
    %eval(sprintf('load %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    eval(sprintf('save %s/%sNSEMfull.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    %close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfullfilters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
    
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfullrasters'], 'Bounds', 'loose', 'Color', 'rgb');
    pause()
    close all

end

%%{
clear
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2016-04-21-1/';
datarun_class = load_data([Analysis_Path 'streamed/data015/data015'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1';
monitor_refresh = 119.5;
cells_orig = [811 3169 2341 6047 4501];
cells_fit = [1053, 3288, 2341, 6049, 4624];
cells_test = [813 3287 2341 6048 4501];
load('/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma2.mat');
mask = imresize(mask,1/4, 'box');
%}

fit_data = 'data019';
test_data = 'data017';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
test_datarun = load_data([Analysis_Path '/mVision/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
stim_length = (3600)*convergence;
fitmovie = uint8(double(fitmovie).*repmat(mask, [1 1 size(fitmovie,3)]) + 64*ones(size(fitmovie)).*(1-repmat(mask, [1 1 size(fitmovie,3)])));
repeats = interleaved_data_prep(test_datarun, 1100, 29, 'cell_spec', cells_test, 'visual_check', 1);
testmovie = fitmovie(:,:,1:1200);
%%
for i = 1
    disp(i)
    glm_cellinfo.cid           = cells_fit(i);
    glm_cellinfo.cell_savename = num2str(cells_fit(i));
    master_idx         = find(fit_datarun.cell_ids == cells_fit(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);    
    fitspikes = fitspikes(fitspikes < stim_length);
    [STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
   
    fittedGLM     = glm_fit(fitspikes, fitmovie,center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
    [gen_signal, firing_rate] = glm_fit_NL_TEST(fittedGLM, fitspikes, fitmovie, 10);
    %[gen_signal, firing_rate] = glm_gen_signal_TEST(fittedGLM, fitmovie);
	model = fitnlm(gen_signal, firing_rate, 'y~b1/(b2+exp(b3*x1))', [1 1 -1]);
    t_init = fittedGLM.linearfilters.Stimulus.time_rk1;
    gen_signal_spatial = glm_gen_signal_spatial(fittedGLM, fitmovie);
    new_timecourse = fminunc(@(time_course)timecourse_error_function(time_course, firing_rate, model.Coefficients.Estimate, gen_signal_spatial), t_init); 

    %eval(sprintf('load %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    eval(sprintf('save %s/%sNSEM.mat fittedGLM', dsave, glm_cellinfo.cell_savename));
    %close all
    plotfilters(fittedGLM);
    set(gcf, 'Position', [100 100 800 250])
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMfilters'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
    
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    exportfig(gcf, [dsave '/' glm_cellinfo.cell_savename '_NSEMrasters'], 'Bounds', 'loose', 'Color', 'rgb');
    close all

end
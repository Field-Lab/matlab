
%%{
clear 
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2016-04-21-1/';
datarun_class = load_data([Analysis_Path 'streamed/data015/data015'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1/NSEM_full_opt';
mkdir(dsave)
monitor_refresh = 119.5;
cells_orig = [811 3169 6047 2281 2341 4501];
cells_fit = [814 3291 6050 2282 2341 4503];
cells_test = [812 3170 6054 2286 2347 4624];
%load('/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma2.mat');
%mask = imresize(mask,1/4, 'box');
%}

fit_data = 'data016';
test_data = 'data018';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
test_datarun = load_data([Analysis_Path '/mVision/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
testmovie = fitmovie(:,:,1:1200);
load('/Volumes/Lab/Users/Nora/downsampledNSbrownian_A_1600.mat');
stim_length = (1500)*convergence;
fitmovie = fitmovie(:,:,1:(1500*120));
repeats = interleaved_data_prep(test_datarun, 1100, 29, 'cell_spec', cells_test, 'visual_check', 1);
datarun_class = load_sta(datarun_class, 'load_sta', cells_orig);

%%
for i = 1:6
    disp(i)
    glm_cellinfo.cid           = cells_fit(i);
    glm_cellinfo.cell_savename = num2str(cells_fit(i));
    master_idx         = find(fit_datarun.cell_ids == cells_fit(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
    fitspikes = fitspikes(fitspikes < stim_length);
    %figure;
    %[STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
    [params,sta,~] = fit_sta(datarun_class.stas.stas{datarun_class.cell_ids == cells_orig(i)}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
    fittedGLM     = glm_fit(fitspikes, fitmovie,round([params.center_point_y, params.center_point_x]), 'monitor_refresh', monitor_refresh, 'WN_STA', squeeze(sum(sta,3)));
    
    %%
    t_init = fittedGLM.linearfilters.Stimulus.time_rk1;
    NL_init = [0.1 0 -4];
    [model, t_init, increments] = fit_NL_and_time(2, fittedGLM, fitspikes, fitmovie, t_init, NL_init);

    % new predictions
    fittedGLM.xvalperformance = glm_predict(fittedGLM, testmovie,'testspikes', repeats.testspikes(:,i));
    new_GS = conv(glm_gen_signal_spatial(fittedGLM, testmovie), t_init, 'full')+fittedGLM.linearfilters.TonicDrive.Filter;
    new_pred = predict(model, new_GS(1:size(testmovie,3))');
    firing =  IDP_plot_PSTH(repeats,1, 'color',0, 'smoothing', 1);
    
    % save the good stuff
    opt_model.orig_glm = fittedGLM;
    opt_model.orig_glm_corr = corr(firing, fittedGLM.xvalperformance.glm_ratepersec(1,1:10:11950)');
    opt_model.raster = Poisson_spiking(new_pred,size(repeats.testspikes,1),fittedGLM.bins_per_frame,monitor_refresh);
    opt_model.corr = corr(firing, new_pred(1:1195));
    opt_model.model = model;
    opt_model.new_time = t_init;
    opt_model.sta_fit.sta = sta;
    opt_model.sta_fit.params = params;
    eval(sprintf('save %s/%sNSEM_optmodel.mat opt_model', dsave, glm_cellinfo.cell_savename));

    %{
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
%}
end
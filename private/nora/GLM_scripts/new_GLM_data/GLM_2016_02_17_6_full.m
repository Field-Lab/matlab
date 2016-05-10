
%%{
clear 
convergence = 1;
Analysis_Path = '/Volumes/Analysis/2016-02-17-6/';
datarun_class = load_data([Analysis_Path 'streamed/data000/data000'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/2016-06-17-6/NSEM_full_opt_refitMU';
mkdir(dsave)
monitor_refresh = 119.5;
cells_orig = [6023 1052 3421  7522 4653  2221 5929 1113 3668 7414 4592 2342];
cells_fit = [6019 1052 3542 7519 4655 2223 5930 1112 3664 7412 4711 2344];
cells_test = [6017;1052;3423;7518;4652;2222;5929;1112;3671;7276;4699;2342];
fit_data = 'data006';
test_data = 'data003';
fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
test_datarun = load_data([Analysis_Path '/mVision/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
load('/Volumes/Lab/Users/Nora/downsampledNSinterval.mat');
testmovie = fitmovie(:,:,1:1200);
fitmovie = fitmovie(:,:,1201:end);
repeats = interleaved_data_prep(test_datarun, 1100, 29, 'cell_spec', cells_test, 'visual_check', 1);
datarun_class = load_sta(datarun_class, 'load_sta', cells_orig);

%%
for i = 1:12
    disp(i)
    glm_cellinfo.cid           = cells_fit(i);
    glm_cellinfo.cell_savename = num2str(cells_orig(i));
    master_idx         = find(fit_datarun.cell_ids == cells_fit(i));
    fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
    fitspikes = fitspikes(fitspikes >(1200/monitor_refresh))-(1200/monitor_refresh);
    fitspikes = fitspikes(fitspikes < (3600-10)*120/monitor_refresh);
    %figure;
    %[STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
    %pause()
    [params,sta,~] = fit_sta(datarun_class.stas.stas{datarun_class.cell_ids == cells_orig(i)}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
    fittedGLM     = glm_fit(fitspikes, fitmovie,round([params.center_point_y, params.center_point_x]), 'monitor_refresh', monitor_refresh, 'WN_STA', squeeze(sum(sta,3)));
    
    %%
    t_init = fittedGLM.linearfilters.Stimulus.time_rk1;
    NL_init = [0.1 0 -4 fittedGLM.linearfilters.TonicDrive.Filter];
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
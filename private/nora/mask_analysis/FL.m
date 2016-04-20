dsave = '/Volumes/Lab/Users/Nora/GLMFits_masking/';
%{
piece = '2016-02-17-1/';
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
%}

%%{
piece = '2015-12-18-2/';
Analysis_Path = '/Volumes/Analysis/2015-12-18-2/data002-data015';
datarun_class = load_data([Analysis_Path '/data002/data002'], struct('load_neurons', 0, 'load_params', 1));
%}

fits_dsave = [dsave piece '/WN/'];

%%
monitor_refresh = 119.5;
cells = [2359 3271 5478 6648]; %167 1082
fit_data = 'data012';
test_datarun = load_data('2015-12-18-2/data002-data015/data013/data013');
test_datarun = load_neurons(test_datarun);
repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 11111);
glm_fit_from_WN(cells, [Analysis_Path '/' fit_data '/' fit_data], 'BW-8-1-0.48-11111', 'testmovie', repeats.testmovie, 'testspikes', repeats.testspikes, 'd_save', fits_dsave, 'monitor_refresh', monitor_refresh, 'stim_length', 1800, 'center_verification', 1);

%%
%test_datarun = load_data('2015-12-18-2/data002-data015/data013/data013');
%test_datarun = load_neurons(test_datarun);
%repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 11111);

%%

%{
n_cells = length(cells);
BPS = zeros(2, n_cells);
for i_cell = 1:n_cells
    disp(i_cell)
    load([dsave '2012-12-18-2/' num2str(cells(i_cell)) '.mat']);
    fittedGLM.xvalperformance = glm_predict(fittedGLM, repeats.testmovie, 'testspikes', repeats.testspikes(:,i_cell));
    plotrasters(fittedGLM.xvalperformance, fittedGLM);
    
    [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis);
    BPS(1,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike/crm_bps;
    load([dsave piece num2str(cells(i_cell)) 'NSEM.mat']);
    [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis);
    BPS(2,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike/crm_bps;
end
figure; plot(BPS(1,:), BPS(2,:), '.')
% save([dsave piece 'BPSOn.mat'], 'BPS')
%}


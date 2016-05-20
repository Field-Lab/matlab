clear
Analysis_Path = '/Volumes/Analysis/2016-01-05-0/';
class = 'data015-mVision/data015/data015';
datarun_class = load_data([Analysis_Path class], struct('load_neurons', 0, 'load_params', 1));
mkdir('/Volumes/Lab/Users/Nora/GLMFits/RPE/201601050/WN/');
cell_type = {'On Parasol'};
cell_type_short = 'OnPar_';

WN_fit = 'mVision/data013/data013';
WN_test = 'mVision/data014/data014';

cells = get_cell_ids(datarun_class, cell_type{1});

%%
ref_data = '/Volumes/Analysis/2016-01-05-0/data015-mVision/data015/';

% fitting 
datarun_fit = load_data([Analysis_Path WN_fit], struct('load_neurons', 1, 'load_params', 1));
new_data = '/Volumes/Analysis/2016-01-05-0/mVision/data013/';
EI_map = crossIdentifyNeuronIDs_NB(ref_data, new_data, cells, [] , 0);
cells_fit = EI_map(:,2);
fitmovie = get_WN_movie('/Volumes/Lab/Users/Nora/Stimulus/BW_XML/BW-16-1-0.48-11111.xml', 120*1800);
monitor_refresh = 100/median(diff(datarun_fit.triggers));

%% testing
datarun_test = load_data([Analysis_Path WN_test], struct('load_neurons', 1, 'load_params', 1));
new_data = '/Volumes/Analysis/2016-01-05-0/mVision/data014/';
EI_map = crossIdentifyNeuronIDs_NB(ref_data, new_data, cells, [] , 0);
cells_test = EI_map(:,2);
prepped_data = interleaved_data_prep(datarun_test, 1100, 30,'cell_spec', cells_test,'visual_check', 0);
testmovie = get_WN_movie('/Volumes/Lab/Users/Nora/Stimulus/BW_XML/BW-16-1-0.48-22222.xml', 1100);

%%
%
for i=1%:length(cells)
    %if sum(cell2mat(prepped_data.testspikes(:,i)))>0.8*10^4
    fit_idx = get_cell_indices(datarun_fit, cells_fit(i));
    fitspikes = align_spikes_triggers(datarun_fit.spikes{fit_idx}, datarun_fit.triggers, 100, monitor_refresh);
    [STA, center] = STA_Test(fitspikes, fitmovie, 1, 1/monitor_refresh);
    fittedGLM = glm_fit(fitspikes, fitmovie, center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
    fittedGLM = save_glm(fittedGLM, prepped_data.testspikes(:,i), testmovie, ['/Volumes/Lab/Users/Nora/GLMFits/RPE/201601050/WN/' cell_type_short num2str(cells(i))]);
    %end
end
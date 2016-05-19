% RPE, WN
clear
Analysis_Path = '/Volumes/Analysis/2014-09-10-1/data004-data013/';
class = 'data012';
datarun_class = load_data([Analysis_Path class '/' class], struct('load_neurons', 0, 'load_params', 1));
mkdir('/Volumes/Lab/Users/Nora/GLMFits/RPE/201409101/WN/');
cell_type = {'On Parasol'};
cell_type_short = 'OnPar_';

WN_run = 'data013';
cells = get_cell_ids(datarun_class, cell_type{1}); % cell ids to fit
cell_idx = get_cell_indices(datarun_class, cell_type{1});
datarun_class = load_sta(datarun_class, 'load_sta', cells);
datarun = load_data([Analysis_Path WN_run '/' WN_run], struct('load_neurons', 1, 'load_params', 1));
prepped_data = interleaved_data_prep(datarun, [3600 1200], 60,'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1','seed', 11111, ...
    'variable_seed_start', 11111);
%prepped_data = interleaved_data_prep(datarun, [3600 1200], 60,'cell_spec', cells,'visual_check', 0);
monitor_refresh = 100/median(diff(datarun.triggers));
fitmovie = concat_movie(prepped_data.fitmovie);

%%
for i=1:length(cells)
    if sum(cell2mat(prepped_data.testspikes(:,i)))>0.8*10^4   
        fitspikes = concat_spikes(prepped_data.fitspikes(:,i), 3600/monitor_refresh);
        [STA, center] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
        fittedGLM = glm_fit(fitspikes, fitmovie, center, 'monitor_refresh', monitor_refresh, 'WN_STA', STA);
        saveGLM(fittedGLM, prepped_data.testspikes(:,i), prepped_data.testmovie, ['/Volumes/Lab/Users/Nora/GLMFits/RPE/201409101/WN/' cell_type_short num2str(cells(i))]);
    end
end

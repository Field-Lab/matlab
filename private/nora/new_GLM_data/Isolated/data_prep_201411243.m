%% isolated 

clear
Analysis_Path = '/Volumes/Analysis/2014-11-24-3/data009-data0012';
datarun_class = load_data([Analysis_Path '/data012/data012'], struct('load_neurons', 0, 'load_params', 1));
cells = get_cell_ids(datarun_class, 'Off Parasol');
%dsave = '/Users/Nora/Desktop/GLMFits/2014-11-24-3';
%mkdir(dsave)
%cells = cells(1);
monitor_refresh = 120;

%% BW
data = 'data011';
datarun = load_data([Analysis_Path '/' data '/' data], struct('load_neurons', 1, 'load_params', 1));
prepped_data = interleaved_data_prep(datarun, [3600 1200], 60, 'cell_spec', cells,'visual_check', 1);%0, 'stimulus_name', 'BW-8-1');

%%
% clear
% Analysis_Path = '/Volumes/Analysis/2014-11-05-6/data002-data005';
% datarun_class = load_data([Analysis_Path '/data002/data002'], struct('load_neurons', 0, 'load_params', 1));
% cells = get_cell_ids(datarun_class, 'On Parasol');
% %dsave = '/Users/Nora/Desktop/GLMFits/2014-11-24-3';
% %mkdir(dsave)
% %cells = cells(1);
% monitor_refresh = 120;
% data = 'data004';
% datarun = load_data([Analysis_Path '/' data '/' data], struct('load_neurons', 1, 'load_params', 1));
% prepped_data = interleaved_data_prep(datarun, [3600 1200], 60, 'cell_spec', cells,'visual_check', 1);%'stimulus_name', 'BW-8-1');

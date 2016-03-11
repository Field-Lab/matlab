clear
Analysis_Path = '/Volumes/Analysis/2014-09-10-2/data002-data005';
datarun_class = load_data([Analysis_Path '/data004/data004'], struct('load_neurons', 0, 'load_params', 1));
cells = get_cell_ids(datarun_class, 'Off Parasol');
dsave = '/Users/Nora/Desktop/GLMFits/2014-09-10-2';
mkdir(dsave)
%cells = cells(1);
monitor_refresh = 120;

%% BW
data = 'data005';
datarun = load_data([Analysis_Path '/' data '/' data], struct('load_neurons', 1, 'load_params', 1));
prepped_data = interleaved_data_prep(datarun, [3600 1200], 60, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1');

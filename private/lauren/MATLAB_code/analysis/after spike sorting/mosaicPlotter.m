
addpath(genpath('~snl-e/matlab-standard/code'))

datarun = load_data('2008-08-27-2/data001-lh/data001-lh');
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));
datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
datarun = get_sta_fits_from_vision(datarun,'all');
figure;plot_rf_summaries(datarun,{1})
figure;plot_rf_summaries(datarun,{3})

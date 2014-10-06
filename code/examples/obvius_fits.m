
% initialize datarun
datarun = load_data('2005-08-08-0','rf-1-auto-jg-0');

% load params file
datarun = load_params(datarun,'verbose',1);

% load neurons file
datarun = load_neurons(datarun);

% load obvius fits from disk
datarun = load_obvius_sta_fits(datarun);

% translate them to matlab format
datarun = get_sta_fits_from_obvius(datarun,'all');

% plot ON parasol fits
figure;plot_rf_summaries(datarun,{1},'plot_fits',1,'label',1)


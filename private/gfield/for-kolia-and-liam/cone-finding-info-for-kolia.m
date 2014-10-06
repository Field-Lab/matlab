% This script loads in data from plantain for Kolia and Liam's project using
% MCMC techniques for Bayesian cone finding

datarun = load_data('2008-08-27-5', 'rf-3-plantain'); % plantian
datarun = load_params(datarun);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', {1,2,3,4,5});

datarun = get_sta_summaries(datarun, {1,2,3,4,5});
datarun = set_polarities(datarun);



%% Analysis for Adaptation Experiment
% authors: Daniel Ahn and Greg Field 
% date: 2011-04-26

% this script loads in data and calculates the mean and standard deviation
% of the receptive field size across midget and parasol RGCs

%% load data

% create data structure: this makes a structure with 'path' information
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-03-25-4/data004/data004/data004');

% load information from the Vision .params file
datarun = load_params(datarun, struct('verbose',1));

% load information from the Vision .sta file
datarun = load_sta(datarun, 'load_sta', 'all');

% load information from the Vision .neurons file
datarun = load_neurons(datarun, 'load_spikes', 'all'); % this isn't necessary, just for illustration

% load information from the Vision .ei file
datarun = load_ei(datarun, 'all', 'array_type', 512);  % this isn't necessary, just for illustration

% we are finished loading data

%% calculate the receptive field size

% Identify cell types of interest
cell_types = {1,2,3,4};

% get sta fits so that we can calculat the rf sizes
datarun = get_sta_fits_from_vision(datarun, cell_types);

% make a list of the RF radii for the on parasol cells
on_parasol_radii = get_rf_fit_radius(datarun, {1});

% calculate the mean of this list
XX

% calculate the standard deviation of this list
XX

% plot a histogram of these values
XX

% repeat the above fot the off parasol, on midget, and off midget cells
XX

% ADVANCED: try to plot the receptive field fits for the on parasol cells.
% HINT: help plot_rf_summaries
XX
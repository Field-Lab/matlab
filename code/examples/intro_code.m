% Intro to the existing lab codes

% first run something like this (to wherever Dropbox is)
% addpath(genpath('Dropbox/Lab/Development/matlabstandard/code'));
% and this one
% Java library
javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

%   Initializing the data structure
%   The most important function is load data! load_data is expecting to be
%   directed to the file with the data000.bin, data000.neuron etc. These are
%   in the Analysis drive. The function is set up so that you can just type
%   'Date/data000'. This initializes the datarun structure. 
datarun=load_data('2012-08-09-3/data002');

% Loading other information
%   Other information you might want including STAs, params, ei, neurons
%   (includes spike times) and polarities
% datarun=load_sta(datarun);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
%datarun=set_polarities(datarun);

% Let's look at the data! This is a classification run.

% Plot the mosaic for a cell type, in this case on parasols
figure;plot_rf_fit(datarun,'On Parasol')
% 
% % Plot the STAs for one cell type, in this case on parasols
% plot_rf_portraits(datarun,{1},'plot_radius',5);

% Look at the firing rate of one cell
get_psth(datarun.spikes{8}, datarun.triggers,'plot_hist', true, 'bin_size',0.05);

% This is a run with rasters.
datarun=load_data('2012-08-09-3/data003');
datarun=load_params(datarun);
datarun=load_neurons(datarun);
% Look at a raster plot of one cell
get_raster(datarun.spikes{8},11*(0:99));

% Take some took to look through the datarun structure to see what all is
% in there. 
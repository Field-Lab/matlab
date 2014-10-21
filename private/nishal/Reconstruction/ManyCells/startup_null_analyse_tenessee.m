%% Startup file for defining matlab and java paths
disp('Running startup file Tennessee...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('~/Dropbox/Lab/Development/matlab/code'));
addpath(genpath('../'));

addpath(genpath('../../plotSpikeRaster_v1'));

addpath(genpath('../../create_act_2/'));

% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')


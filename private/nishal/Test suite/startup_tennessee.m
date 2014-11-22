%% Startup file for defining matlab and java paths
disp('Running startup file ');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/Users/vision/GITs/matlab/code'));
addpath(genpath('/Users/vision/GITs/matlab/private/nora'));

addpath(genpath('/Users/vision/GITs/matlab/private/nishal/plotSpikeRaster_v1'));
addpath(genpath('/Users/vision/GITs/matlab/private/nishal/create_act2'));

% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

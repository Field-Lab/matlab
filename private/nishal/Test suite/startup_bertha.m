%% Startup file for defining matlab and java paths
disp('Running startup file ');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/home/vision/Nishal/matlab/code'));
addpath(genpath('/home/vision/Nishal/matlab/private/nora'));

addpath(genpath('/home/vision/Nishal/matlab/private/nishal/plotSpikeRaster_v1'));
addpath(genpath('/home/vision/Nishal/matlab/private/nishal/create_act_2'));
% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')
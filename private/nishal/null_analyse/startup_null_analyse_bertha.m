%% Startup file for defining matlab and java paths
disp('Running startup file Bertha...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('~/Nishal/matlab/code'));

% set some default plot 
set(0, 'DefaultAxesFontSize', 14, 'DefaultAxesFontName', 'Helvetica')


addpath(genpath('/Volumes/Lab/Users/bhaishahster/cvx'))
cvx_setup
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/plotSpikeRaster_v1'));

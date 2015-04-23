%% Startup file for defining matlab and java paths
disp('Running startup file Bertha...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('~/Nishal/matlab/code'));





addpath(genpath('~/Nishal/matlab/private/nishal/plotSpikeRaster_v1'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/'))
% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')


addpath(genpath('/Volumes/Lab/Users/bhaishahster/cvx'))
cvx_setup
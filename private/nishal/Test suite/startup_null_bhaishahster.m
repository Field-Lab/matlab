%% Startup file for defining matlab and java paths
disp('Running startup file ');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Users/bhaishahster/GITs/matlab/private/nora'));
addpath(genpath('/Users/bhaishahster/GITs/matlab/private/lauren/MATLAB_code')); 
addpath(genpath('/Users/bhaishahster/GITs/private/nishal/plotSpikeRaster_v1'));

% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

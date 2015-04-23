%% Startup file for defining matlab and java paths
disp('Running startup file Tennessee...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('~/GITs/matlab/code'));

addpath(genpath('~/GITs/matlab/private/lauren/MATLAB_code')); 



addpath(genpath('~/GITs/matlab/private/nishal/plotSpikeRaster_v1'));

addpath(genpath('~/GITs/matlab/private/nishal/Test suite'));
% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')


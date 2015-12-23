%% Startup file for defining matlab and java paths
disp('Running startup file Tennessee...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/code'));
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/private/nora'));
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/private/nora/generalcomputations_AKH'));

addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/private/lauren/MATLAB_code')); 

% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

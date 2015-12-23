%% Startup file for defining matlab and java paths
disp('Running startup file 5/27/14...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/code'));
addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/private/nora'));
addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/private/nora/generalcomputations_AKH'));

addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/private/lauren/MATLAB_code')); 

% set some default plot 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

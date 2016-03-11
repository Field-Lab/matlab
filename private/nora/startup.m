%% Startup file for defining matlab and java paths
disp('Running startup file 1/27/14...');

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths
addpath(genpath('/home/vision/Nora/matlab/code/projects/glm'));
addpath(genpath('/home/vision/Nora/matlab/code/lab'));
addpath(genpath('/home/vision/Nora/matlab/private/nora'));


% set some default plot stuff
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')



% addpath(genpath('/home/vision/Nishal/matlab/code'));

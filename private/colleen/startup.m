% Startup file for defining matlab and java paths
disp('Running startup file...');

% Java library
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths 
addpath(genpath('/Users/vision/Desktop/GitHub code repository/code'));
% addpath(genpath('~/matlab'));

% set some default plot stuff
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

fprintf('Running startup file from %s\n',mfilename('fullpath'));
% Java library
javaaddpath('/Volumes/Lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths 
addpath(genpath('/Volumes/Lab/Users/grosberg/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/grosberg/matlab/private/lgrosberg'));
addpath(genpath('/Volumes/Lab/Users/grosberg/direwolf_matlab')); 
% cd /Volumes/Lab/Users/grosberg/cvx
% cvx_setup /Volumes/Lab/Users/grosberg/cvx_license.dat
run /Volumes/Lab/Users/grosberg/cvx/cvx_startup.m

set(0,'DefaultAxesFontSize',18,'DefaultAxesFontName','Myriad Pro');
set(0,'DefaultFigureColor','w'); 
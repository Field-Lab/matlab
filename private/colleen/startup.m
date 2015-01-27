% Startup file for defining matlab and java paths
disp('Running startup file...');

addpath('/home/vision/alex/matlab/code/lab/utilities')
addpath(genpath2('/home/vision/alex/matlab/code', {'.svn'}));
visstable()
%addpath('/home/vision/alex/matlab/')
%addpath(genpath('/home/vision/alex/matlab/private/alex/retina-subunits-master'))

% Java library
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths 
addpath(genpath('/Users/colleen/matlab/code'));
addpath(genpath('~/matlab'));

% Add utilities so that you have genpath2.m available
addpath(genpath('/Users/colleen/matlab/utilities'))

% Use genpath2 to recursively add everything in matlab-standard to your path, but exclude .svn directories
addpath(genpath(['../../../code']));

% Local Java builds: Vision, Cell-Finder
% We'll have to fix vision_path_stable function to point to the right place on wharf-d
visstable();

% set some default plot stuff
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')


addpath(genpath('/Users/colleen/matlab/private/colleen'));



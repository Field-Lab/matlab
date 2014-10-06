disp('Running startup...');
% Java library
javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');

% Add utilities so that you have genpath2.m available
addpath([ '../../../code/lab/utilities'])

% Use genpath2 to recursively add everything in matlab-standard to your path, but exclude .svn directories
addpath(genpath2(['../../../code'], {'.svn'}));

addpath(genpath(['../cvx']))
cvx_setup

addpath(genpath(['../CBPSpikesortDemoPackage']))

% Local Java builds: Vision, Cell-Finder
% We'll have to fix vision_path_stable function to point to the right place on wharf-d
visstable();
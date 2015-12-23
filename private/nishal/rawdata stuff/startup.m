disp('Running startup...');

% Add utilities so that you have genpath2.m available
addpath([ '../../../code/lab/utilities'])

% Use genpath2 to recursively add everything in matlab-standard to your path, but exclude .svn directories
addpath(genpath2(['../../../code'], {'.svn'}));

% Local Java builds: Vision, Cell-Finder
% We'll have to fix vision_path_stable function to point to the right place on wharf-d
visstable();
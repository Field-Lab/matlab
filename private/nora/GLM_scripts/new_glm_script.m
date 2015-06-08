d% addpath(genpath('/Users/Nora/Documents/MATLAB/matlab/code/projects/glm'))
addpath(genpath('/home/vision/Nora/matlab/code/projects/glm')) %BERTHA

glmpath_Lovelight

%% Common changes

% changes{1}.type = 'filter_mode';
% changes{1}.name ='rk2';

% changes{1}.type = 'debug';
% changes{1}.name ='true';

changes{1}.type = 'CouplingFilters';
changes{1}.name ='ON';

%% Other inputs
experiments = [1 3 4]; % 1-4
stimulus = [1 2]; % 1 is WN, 2 is NSEM
celltypes = [2]; % 1 is ON, 2 is OFF
cellsubset = 'all'; % options are debug, shortlist, or all;

%% Run glmwrap
glm_wrap(experiments,stimulus,celltypes,cellsubset,changes)
glm_wrap(2,1,celltypes,cellsubset,changes)

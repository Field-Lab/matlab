addpath(genpath('/home/vision/Nora/matlab/code/projects/glm')) %BERTHA

glmpath_Lovelight


%% Set up the different fit types I want to test
% Let's just start with rk1 and coupling separately.

% Basic fit convergence
% leaving changes{1} blank

% Rk1 convergence
changes{2}{1}.type = 'filter_mode';
changes{2}{1}.name ='rk1';

% Rk2 convergence
changes{3}{1}.type = 'filter_mode';
changes{3}{1}.name ='rk2';

%% Which cells and stimuli to test

experiments = [1 2 3 4]; % 1-4
stimulus = [1 2]; % 1 is WN, 2 is NSEM
celltypes = [1 2]; % 1 is ON, 2 is OFF
cellsubset = 'debug'; % options are debug, shortlist, or all;

%% Run glmwrap
for i = 1:3
    glm_wrap_convergence(experiments,stimulus,celltypes,cellsubset,changes{i})
end
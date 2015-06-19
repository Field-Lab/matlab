addpath(genpath('/home/vision/Nora/matlab/code/projects/glm')) %BERTHA

glmpath_Lovelight


%% Set up the different fit types I want to test
% Let's just start with rk1 and coupling separately.

% Basic fit convergence
% leaving changes{1} blank

% % coupling convergence
% changes{2}.type = 'CouplingFilters';
% changes{2}.name = 'ON';

% cone model
changes{2}.type = 'cone_model';
changes{2}.name = 'rieke_linear';

% 
% % Rk1 convergence
% changes{2}{1}.type = 'filter_mode';
% changes{2}{1}.name ='rk2';

% Rk2 convergence
changes{1}.type = 'filter_mode';
changes{1}.name ='rk2';

%% Which cells and stimuli to test

experiments = [1 2 ]; % 1-4
stimulus = [2]; % 1 is WN, 2 is NSEM
celltypes = [1 2]; % 1 is ON, 2 is OFF
cellsubset = 'debug'; % options are debug, shortlist, or all;

%% Run glmwrap
for i = 1
    glm_wrap_convergence(experiments,stimulus,celltypes,cellsubset,changes)
end

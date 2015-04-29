addpath(genpath('/home/vision/Nora/matlab/code/projects/glm')) %BERTHA
glmpath_Lovelight
changes{1}.type = 'filter_mode';
changes{1}.name = 'rk1';
glm_init(1, 1, [1 2], 1, 'debug', changes);
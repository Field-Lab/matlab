% addpath(genpath('/Users/Nora/Documents/MATLAB/matlab/code/projects/glm'))
addpath(genpath('/home/vision/Nora/matlab/code/projects/glm')) %BERTHA
glmpath_Lovelight

%% Common changes

changes{1}{1}.type = 'filter_mode';
changes{1}{1}.name ='rk1';

changes{1}{2}.type = 'Subunits';
changes{1}{2}.name ='ON';

changes{1}{3}.type = 'debug';
changes{1}{3}.name ='false';

changes{2}{1}.type = 'CouplingFilters';
changes{2}{1}.name ='ON';

changes{2}{1}.type = 'CouplingFilters';
changes{2}{1}.name ='OFF';

changes{3}{2}.type = 'filter_mode';
changes{3}{2}.name ='rk2';

changes{4}{2}.type = 'filter_mode';
changes{4}{2}.name ='rk2';

changes{3}{1}.type = 'CouplingFilters';
changes{3}{1}.name ='ON';

changes{4}{1}.type = 'CouplingFilters';
changes{4}{1}.name ='OFF';

runoptions.replace_existing = true;


%% Other inputs
experiments = [1]; % 1-4
stimulus = [2]; % 1 is WN, 2 is NSEM
celltypes = [1]; % 1 is ON, 2 is OFF
cellsubset = 'all'; % options are debug, shortlist, or all;


for changes_index = 1
     glm_wrap(1,1,1,'all', changes{changes_index}, runoptions);
end

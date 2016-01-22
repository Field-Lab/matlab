% NBrackbill 2015-04-20
% Calls NSEM_BaseDirectories
%%
function glmpath_Lovelight

copyfile('/home/vision/Nora/matlab/code/projects/glm/NSEM_BaseDirectories_Lovelight.m','/home/vision/Nora/matlab/code/projects/glm/NSEM_BaseDirectories.m') %BERTHA
%copyfile('/Users/Nora/Documents/MATLAB/matlab/code/projects/glm/NSEM_BaseDirectories_Lovelight.m' , '/Users/Nora/Documents/MATLAB/matlab/code/projects/glm/NSEM_BaseDirectories.m')
BD = NSEM_BaseDirectories;
glmwrapdir     = BD.GLM_codehome ;
glmcodedir     = sprintf('%s/glm_core', glmwrapdir);
%%
path(sprintf('%s', glmwrapdir), path)
path(sprintf('%s/wrap_bookkeep', glmwrapdir), path)
path(sprintf('%s', glmcodedir), path)
% OLDER VERSIONS REQUIRED JAVA PATH
%run('/home/vision/alex/matlab/private/alex/startup.m')
end



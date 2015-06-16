% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories
%%
function glmpath_Bertha

copyfile('NSEM_BaseDirectories_Bertha.m' , 'NSEM_BaseDirectories.m')
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



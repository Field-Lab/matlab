% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories

function glmpath_Alligator
clear; clc;
%corelabpath_string = alligator_corelab_path;
%addpath(corelabpath_string)
%javaaddpath('/Users/akheitman/Dropbox/Lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
%%
copyfile('NSEM_BaseDirectories_Alligator.m' , 'NSEM_BaseDirectories.m')
BD = NSEM_BaseDirectories;
glmwrapdir     = BD.GLM_codehome ;
glmcodedir     = sprintf('%s/glm_core', glmwrapdir);
%%

path(sprintf('%s', glmwrapdir), path)
path(sprintf('%s/wrap_bookkeep', glmwrapdir), path)
path(sprintf('%s', glmcodedir), path)

end


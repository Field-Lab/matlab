% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories

clear all; close all;

restoredefaultpath;
corelabpath_string = corelab_path;
addpath(corelabpath_string)

javaaddpath('/Users/akheitman/Dropbox/Lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
%%
copyfile('/Users/akheitman/Matlab_code/glm_AH/NSEM_BaseDirectories_Alligator.m' , '/Users/akheitman/Matlab_code/glm_AH/NSEM_BaseDirectories.m')
BD = NSEM_BaseDirectories;
glmcodedir     = BD.GLM_codehome ;
generalcodedir = BD.general_codehome;
%%
% silly path stuff
path(sprintf('%s', glmcodedir), path)
path(sprintf('%s/RunGLM', glmcodedir), path)
path(sprintf('%s/prep', glmcodedir), path)
path(sprintf('%s/prep/conemodels', glmcodedir), path)
path(sprintf('%s/bookkeep', glmcodedir), path)
path(sprintf('%s/newcorecode', glmcodedir), path)
path(sprintf('%s/corecode', glmcodedir), path)
path(sprintf('%s/newtestcode', glmcodedir), path)
path(sprintf('%s/Metrics_Evaluation', glmcodedir), path)
path(sprintf('%s/chaitu/mfiles', glmcodedir), path)
path(sprintf('%s/cellselection', glmcodedir), path)
path(sprintf('%s/generalcomputations_AKH', generalcodedir), path)
path(sprintf('%s/External_Code', generalcodedir), path)
%%
clear all; close all;

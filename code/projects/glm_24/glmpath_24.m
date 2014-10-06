% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories

clear all; close all;

restoredefaultpath;
corelabpath_string = corelab_path;
addpath(corelabpath_string)



BD = NSEM_BaseDirectories;

glmcodedir     = BD.GLM_codehome ;
generalcodedir = BD.general_codehome;
% silly path stuff
path(sprintf('%s', glmcodedir), path)
path(sprintf('%s/RunGLM', glmcodedir), path)
path(sprintf('%s/prep', glmcodedir), path)
path(sprintf('%s/prep/conemodels', glmcodedir), path)
path(sprintf('%s/bookkeep', glmcodedir), path)
path(sprintf('%s/newcorecode', glmcodedir), path)
path(sprintf('%s/corecode', glmcodedir), path)
path(sprintf('%s/Metrics_Evaluation', glmcodedir), path)
path(sprintf('%s/chaitu/mfiles', glmcodedir), path)
path(sprintf('%s/generalcomputations_AKH', generalcodedir), path)


clear all; close all;

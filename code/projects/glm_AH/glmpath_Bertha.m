% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories
%%
clear all; close all;


addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/code'))
labdcodedropbox ='/home/vision/Dropbox/Lab/Development/matlab-standard/code';

% add java path   2015-01-26
display('adding java path from startup file in "/home/vision/alex/matlab/private/alex/startup.m"') 
run('/home/vision/alex/matlab/private/alex/startup.m')
%%

%
%restoredefaultpath;
%corelabpath_string = corelab_path_Bertha;
%addpath(corelabpath_string)
%}

copyfile('/home/vision/akheitman/matlab/code/glm_AH/NSEM_BaseDirectories_Bertha.m' , '/home/vision/akheitman/matlab/code/glm_AH/NSEM_BaseDirectories.m')
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



clear all; close all;

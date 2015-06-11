% AK Heitman 2014-04-02
% Calls NSEM_BaseDirectories

clear all; close all;

restoredefaultpath;
% corelabpath_string = corelab_path;
% addpath(corelabpath_string)
% NB 06-11-2014

% BD = NSEM_BaseDirectories;

glmcodedir = '/Users/colleen/matlab/private/nora/glm_24/';
generalcodedir = '/Users/colleen/matlab/private/nora/';

% glmcodedir     = BD.GLM_codehome ;
% generalcodedir = BD.general_codehome;
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

addpath('../../utilities')

% addpath(genpath2('/home/vision/alex/matlab/code', {'.svn'}));
% addpath(genpath2('matlab/private/alex/matlab/code', {'.svn'}));

% visstable()
%addpath('/home/vision/Colleen/matlab/utilities')
%addpath(genpath2('/home/vision/Colleen/matlab/code', {'.svn'}));
%visstable()
%addpath('/home/vision/alex/matlab/')
%addpath(genpath('/home/vision/alex/matlab/private/alex/retina-subunits-master'))

% Java library
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths 
addpath(genpath('/Users/colleen/matlab/code'));
addpath(genpath('~/matlab'));

% Add utilities so that you have genpath2.m available
addpath(genpath('/Users/colleen/matlab/utilities'))

% Use genpath2 to recursively add everything in matlab-standard to your path, but exclude .svn directories
addpath(genpath(['../../../code']));

% Local Java builds: Vision, Cell-Finder
% We'll have to fix vision_path_stable function to point to the right place on wharf-d
visstable();

% set some default plot stuff
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')
% set figure background to white
get(0, 'Factory');
set(0, 'defaultfigurecolor', [1,1,1]);

addpath(genpath('/Users/colleen/matlab/private/colleen'));

clear all; close all;

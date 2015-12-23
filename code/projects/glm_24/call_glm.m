% Hack storage of command calls to run GLM in Bertha

clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS


% SETUP cells and experiments, the TYPE of GLM (GLMType) 

GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'WN'; GLMType.map_type = 'mapPRJ';
GLMType.debug = true; 
GLMType.specialchange = false;
%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
%GLMType.input_pt_nonlinearity      = true;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.CONVEX = true;
%GLMType.DoubleOpt = true;
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}
GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
BD = NSEM_BaseDirectories;
exptests = [1];
cellselectiontype = 'debug'; %cellselectiontype = 'shortlist';
GLMType.fitname  = GLM_fitname(GLMType);   
glmwrap24_func(GLMType,BD, exptests,cellselectiontype)



%%
clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS


% SETUP cells and experiments, the TYPE of GLM (GLMType) 

GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = true; 
GLMType.specialchange = false;
%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';

GLMType.CONVEX = true;
%GLMType.DoubleOpt = false
GLMType.DoubleOpt = true;
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}
GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
BD = NSEM_BaseDirectories_Bertha;
exptests = [1];
%cellselectiontype = 'debug'; 
cellselectiontype = 'shortlist';
GLMType.fitname  = GLM_fitname(GLMType);   
glmwrap24_func(GLMType,BD, exptests,cellselectiontype)
%%
clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS


% SETUP cells and experiments, the TYPE of GLM (GLMType) 

GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.specialchange = false;
%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';

GLMType.CONVEX = true;
%GLMType.DoubleOpt = false
GLMType.DoubleOpt = true;
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}
GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
BD = NSEM_BaseDirectories_Bertha;
exptests = [1];
cellselectiontype = 'debug'; 
%cellselectiontype = 'shortlist';
GLMType.fitname  = GLM_fitname(GLMType);   
glmwrap24_func(GLMType,BD, exptests,cellselectiontype)
%%
clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
% SETUP cells and experiments, the TYPE of GLM (GLMType) 


GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname ='p8Mod1Max1e4p8';
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.nullpoint = 'mean'; 
GLMType.fit_type = 'WN'; GLMType.map_type = 'mapPRJ';
GLMType.debug = true; 
GLMType.specialchange = false;
%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';

GLMType.CONVEX = true;
GLMType.DoubleOpt = true;
%GLMType.DoubleOpt = true;
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
BD = NSEM_BaseDirectories_Bertha;
exptests = [1];
cellselectiontype = 'debug'; 
%cellselectiontype = 'shortlist';
GLMType.fitname  = GLM_fitname(GLMType);   
glmwrap24_func(GLMType,BD, exptests,cellselectiontype)






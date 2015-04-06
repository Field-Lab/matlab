% AKHeitman 2014-03-27
% AKHeitman 2014-04-07  -emphasize role as parameter independet loading
%                       -setting directories
%                       -dicatating parameters
%                       -no computations
% Calls:
% GLM_fitname
% cell_list
% Directories_Params_v23
% NSEM_secondaryDirectories
% loadmoviematfile
% concat_fitmovie_fromblockedcell
% findcelltype
% concat_fitspikes_fromorganizedspikes
% visionSTA_to_xymviCoord
%%%
% glm_execute
%%%%

% Fit Dependent
% Parameter Independent
% Organize fit movie and spike times, directories etc.
% DO ALL LOADING AND SAVING HERE!!!
% All stimulus parameters and blocks everything should get dealt with here


% movie, spikes, block structure, GLMType
% the Convtest calling sequence needs to get worked out.. but otherwise ok
clear; close all;  clc

%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE!

% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS

% SETUP cells and experiments, the TYPE of GLM (GLMType)

GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean';
GLMType.nullpoint = 'mean';
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false;
GLMType.specialchange = false;
GLMType.specialchange_name = 'extra_coupling';
GLMType.CBP=false;

GLMType.stimfilter_mode = 'rk1';
%GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = false;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';

GLMType.CONVEX = true;
GLMType.DoubleOpt = false;
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
GLMType.Subunits = false;
GLMType.Saccades = false;
GLMType.color = false;

% GLMType.fixed_spatialfilter = true;
% NBCoupling 06-12-2014
GLMType.func_sname = 'glmwrap24_CP';
GLMType.fullmfilename =mfilename('fullpath');
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);
troubleshoot.doit    = false;
%troubleshoot.plotdir = '/Users/akheitman/Matlab_code/troubleshooting_plots'
% NBCoupling
%troubleshoot.plotdir=BD.GLM_troubleshootplots;
troubleshoot.name    = 'singleopt';

%  LOOP THROUGH DATA SETS

BD = NSEM_BaseDirectories;

exptests = [1 2 3 4];
cellselectiontype = 'debug';
troubleshoot.plotdir = BD.GLM_troubleshootplots;
%%
count = 0;
BPS_rk1= zeros(7,2);
for i_exp = exptests
    %%
    expnumber = i_exp;
    [exp_nm,cells,~]  = cell_list( expnumber, cellselectiontype);
    [StimulusPars, ~, ~, datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    
    inputs.exp_nm    = exp_nm;
    inputs.map_type  = GLMType.map_type;
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
    
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs;
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));
    display(sprintf('Save Directory :  %s', d_save));
    %     if ~exist(d_save), mkdir(d_save); end
    GLMType.d_save = d_save;
    
    %% Load Movie and Concatenate the Fitting Section
    [testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');    
    
    for i_cell = 1:length(cells)
        clear glm_cellstruct
        cid = cells{i_cell};
        [~ , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);

        % Load the correct GLM
        try
            load([d_save '/' cell_savename '.mat'])
            disp('Fitted GLM Loaded');
            count = count+1;
            BPS_rk1(count,1) = fittedGLM.xvalperformance.glm_normedbits;
            BPS_rk1(count,2) = cid ;
%             params = fittedGLM.rawfit.opt_params;
%             t1_idx = fittedGLM.rawfit.paramind.time1;
%             t2_idx = fittedGLM.rawfit.paramind.time2;
%             s1_idx = fittedGLM.rawfit.paramind.space1;
%             s2_idx = fittedGLM.rawfit.paramind.space2;
%             time = (1:length(t1_idx))*1/120*1000;
%             figure('Position', [100 100 700 200]);
%             subplot(1,3,1)
%             plot(time, params(t1_idx))
%             hold on
%             plot(time, params(t2_idx))
%             title(cell_savename)
%             xlim([0 250])
%             xlabel('Time (ms)')
%             subplot(1,3,2)
%             imagesc(reshape(params(s1_idx),11,11))
%             caxis([0 0.5])
%             axis image
%             subplot(1,3,3)
%             imagesc(reshape(params(s2_idx),11,11))
%             caxis([0 0.5])
%             axis image
%             print(['Users/Nora/Desktop/' cell_savename '.eps'], 'eps')
        catch
        end
    end
    
end


%%


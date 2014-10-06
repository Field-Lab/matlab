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



%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE!
clear; close all;  clc
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

%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = true;
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
GLMType.CONVEX = true;
GLMType.DoubleOpt = true;
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}

GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = true;
GLMType.fixed_spatialfilter = true;
% NBCoupling 06-12-2014
GLMType.func_sname = 'glmwrap_24_CP';
GLMType.fullmfilename =mfilename('fullpath');
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);
troubleshoot.doit    = true;
%troubleshoot.plotdir = '/Users/akheitman/Matlab_code/troubleshooting_plots'
% NBCoupling
%troubleshoot.plotdir=BD.GLM_troubleshootplots;
troubleshoot.name    = 'singleopt';

%  LOOP THROUGH DATA SETS

BD = NSEM_BaseDirectories;
cellselectiontype = 'debug';
troubleshoot.plotdir = BD.GLM_troubleshootplots

expnumber = 1;
[exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype)
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
cid = cells{2};

% NBCoupling 06-12-14
if GLMType.CouplingFilters==true
    CTYPE = {'On-Parasol','Off-Parasol'};
    for n_ctype = 1:length(CTYPE)
        [~, parasol_cp.ids{n_ctype}] = celltype_id_AH(CTYPE{n_ctype}, datarun_slv.cell_types);
    end
    parasol_cp.indices{1}=get_cell_indices(datarun_mas,parasol_cp.ids{1});
    parasol_cp.indices{2}=get_cell_indices(datarun_mas,parasol_cp.ids{2});
    parasol_cp.NumCells{1}=length(parasol_cp.ids{1});
    parasol_cp.NumCells{2}=length(parasol_cp.ids{2});
end
% end NBCoupling

%%%%  Shorten Block count if using Debug
if GLMType.debug
    StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:10);
end
clear boolean_debug map_type fit_type shead_cellID expname

%%%%%% Name and Create a Save Directory %%%%%%%%%%%

inputs.exp_nm    = exp_nm;
inputs.map_type  = GLMType.map_type;
inputs.stim_type = GLMType.fit_type;
inputs.fitname   = GLMType.fitname;


%% Load Cell Specific Elements   Spikes and STA
inputs.exp_nm       = exp_nm;
inputs.map_type     = GLMType.map_type;
DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs);
inputs.stim_type    = GLMType.fit_type;
DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
clear inputs




[celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);

glm_cellinfo.cid           = cid;
glm_cellinfo.exp_nm        = exp_nm;
glm_cellinfo.celltype      = celltype;
glm_cellinfo.cell_savename = cell_savename;
glm_cellinfo.fitname       = GLMType.fitname
glm_cellinfo.computedtstim = StimulusPars.slv.computedtstim;
%%% only really matters when binning spikes for fitting %%
%%% small differences from 1/120 will start adding up over the long run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load the center_coord  (do Param dependent ROI manipulation later
master_idx         = find(datarun_mas.cell_ids == cid);
stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
slvdim.height      = StimulusPars.slv.height; slvdim.width = StimulusPars.slv.width;
[center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
glm_cellinfo.slave_centercoord = center_coord;
glm_cellinfo.pairs=pick_neighbor_cells_test(stafit_centercoord, parasol_cp, datarun_mas.vision.sta_fits);

figure;
title('Coupled Cells')
hold on
plot_rf_fit(datarun_mas, parasol_cp.ids{1},'edge_color',[0 1 0])
plot_rf_fit(datarun_mas, parasol_cp.ids{2},'edge_color',[1 0 0])
plot_rf_fit(datarun_mas, glm_cellinfo.pairs(1:6),'fill_color',[0 1 0],'fill',true,'edge',false)
plot_rf_fit(datarun_mas, glm_cellinfo.pairs(7:12), 'fill_color',[1 0 0],'fill',true,'edge',false)
plot_rf_fit(datarun_mas, glm_cellinfo.cid,'fill',true)

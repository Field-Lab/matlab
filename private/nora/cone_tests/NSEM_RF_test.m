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
GLMType.CBP=false;

%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = false;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';

GLMType.CONVEX = true;
GLMType.DoubleOpt = false;

GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.Subunits = false;
GLMType.Saccades=false;
GLMType.color=false;

% GLMType.fixed_spatialfilter = true;
% NBCoupling 06-12-2014
GLMType.func_sname = 'glmwrap24_CP';
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

exptests = [3];
cellselectiontype = 'shortlist';
troubleshoot.plotdir = BD.GLM_troubleshootplots
%%

%%
expnumber = i_exp;
[exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);

clear boolean_debug map_type fit_type shead_cellID expname

inputs.exp_nm    = exp_nm;
inputs.map_type  = GLMType.map_type;
inputs.stim_type = GLMType.fit_type;
inputs.fitname   = GLMType.fitname;

%% Load Movie and Concatenate the Fitting Section
clear Main_SolPars Other_SolParss
%%% Load Stimulus   -- insert more frame cutting here!
[blockedmoviecell, ~, ~] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
[testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');


%% Load Cell Specific Elements   Spikes and STA
inputs.exp_nm       = exp_nm;
inputs.map_type     = GLMType.map_type;
DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs);
inputs.stim_type    = GLMType.fit_type;

DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
clear inputs

clear glm_cellstruct
cid = cells{1};
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
eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
spikesconcat.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
eval(sprintf('load %s/STAandROI_%s.mat STAandROI', DirPars.WN_STAdir, cell_savename));
glm_cellinfo.WN_STA = STAandROI.STA;
% clear cell_savename BD DirPars GLMType StimulusPars cellselectiontype celltype datarun_mas datarun_slv exp_nm expnumber exptests glm_cellinfo

%% Calculate the STA
STA=zeros(80,40,30);
stim = double(testmovie{1}.matrix);
test_STA=1;
for i_block=1:20
    spikes = organizedspikes.block.t_sp{2*i_block-test_STA};
    stim = double(blockedmoviecell{i_block}.matrix);
    frame_times = organizedspikes.block.t_frame{2*i_block-test_STA};
    for i_spike=1:length(spikes)
        [~,spike_frame]=min(abs(frame_times-spikes(i_spike)));
        if spike_frame>29
            STA=STA+stim(:,:,spike_frame-29:spike_frame);
        end
    end
end
%%

for i=1:30
    subplot(2,1,1)
    imagesc(squeeze(STA(:,:,i))');
    axis image
    colormap gray
    hold on; scatter(center_coord.x_coord, center_coord.y_coord); hold off;
    title('NSEM STA')
    
    subplot(2,1,2)
    imagesc(squeeze(STAandROI.STA(:,:,31-i))');
    axis image
    colormap gray
    hold on; scatter(center_coord.x_coord, center_coord.y_coord); hold off;
    title('White Noise STA')
    
    pause()
end

%%
i=27;
subplot(2,1,1)
imagesc(squeeze(STA(:,:,i))');
axis image
colormap gray
hold on; scatter(center_coord.x_coord, center_coord.y_coord); hold off;
title('NSEM STA')

subplot(2,1,2)
imagesc(squeeze(STAandROI.STA(:,:,31-i))');
axis image
colormap gray
hold on; scatter(center_coord.x_coord, center_coord.y_coord); hold off;
title('White Noise STA')

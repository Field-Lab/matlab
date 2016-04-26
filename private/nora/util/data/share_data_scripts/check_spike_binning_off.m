% NB 2014-10-9
% get data into a reasonable shape to share

% Data Type
clear; close all;  clc
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.map_type = 'mapPRJ';
i_exp = 2; i_cell = 1;
cellselectiontype = 'shortlist';
GLMType.fit_type = 'NSEM';
tstim = .00832750;

% load basic info
expnumber = i_exp;
[exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
cells
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
clear boolean_debug map_type fit_type shead_cellID expname
inputs.exp_nm    = exp_nm;
inputs.map_type  = GLMType.map_type;
inputs.stim_type = GLMType.fit_type;

% load movie
clear Main_SolPars Other_SolParss
[blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
[testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');
clear origmatfile

for i=1:length(blockedmoviecell)
    NSEMStimData.FitMovie{i}=blockedmoviecell{i}.matrix;
end
NSEMStimData.TestMovie=testmovie{1}.matrix;
clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
%}

% Load STA and spikes
inputs.exp_nm       = exp_nm;
inputs.map_type     = GLMType.map_type;
DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs);
inputs.stim_type    = GLMType.fit_type;

DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
clear inputs
fit_dir = '/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/rk1_MU_PS_noCP_timekernelCONEMODEL_stimnonlin_log_powerraise/standardparams/PS_netinhibitory_domainconstrain_COB/postfilterNL_Logistic_2Par_fixMU/NSEM_mapPRJ/';

%%

for i_cell = 2%:length(cells)
    cid = cells{i_cell};
    %[~ , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);
    cell_savename = 'OFFPar_31';
    master_idx         = find(datarun_mas.cell_ids == cid);
    stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
    stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
    slvdim.height      = StimulusPars.slv.height; slvdim.width = StimulusPars.slv.width;
    [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
    clear master_idx stafit_centercoord slvdim sd
    
    %eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
    load(['/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_' cell_savename '.mat']);
    logical_spike = spike_bin_test(organizedspikes, StimulusPars, NSEMStimData.TestMovie);
    load([fit_dir exp_nm '/' cell_savename '.mat'])

    imagesc(logical_spike); 
    figure; imagesc(fittedGLM.xvalperformance.rasters.recorded);
end


%% loading and organizing stimulus
load(['/Volumes/Lab/Users/Nora/data_files/Data/CarlosData/NSEM-' exp_nm '-CellData.mat'])
eval(['ella_spikes.block.t_sp_withinblock = NSEMCellData.' cell_savename '.Spikes;'])
logical_spike_ella = spike_bin_test(ella_spikes, StimulusPars, NSEMStimData.TestMovie);



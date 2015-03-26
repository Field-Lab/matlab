%%


clear; close all; clc
BD = NSEM_BaseDirectories;
basefolder   = sprintf('%s/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots', BD.NSEM_home)
shortlist = true;
if exist('shortlist','var') && shortlist
    eval(sprintf('load %s/Raw_Conv_shortlist.mat', basefolder));
else
    eval(sprintf('load %s/Raw_Conv.mat', basefolder)); 
end
pct_change = [5:10:95];pct_change = [pct_change, 100];
celltypes = 1:2;
exps = 1:4;
i_exp = 1; i_cell = 1; i_celltype = 1; i_pct = 1;  RCD = Raw_Conv_Data;

debug = true;
%%

% CYCLE THROUGH CELLS 
for i_exp = exps
    %%
    % LOAD DIRECTORIES OF STA AND ORGANIZED SPIKES
    exp_nm  = RCD{i_exp}.exp_nm;
    GLMType = RCD{i_exp}.GLMType;    
    secondDir.exp_nm    = exp_nm; 
	secondDir.map_type  = GLMType.map_type; 
	secondDir.stim_type = 'WN';
	Dirs.WN_organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
	secondDir.stim_type = 'NSEM';
	Dirs.NSEM_organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
    Dirs.WN_STAdir               = NSEM_secondaryDirectories('WN_STA', secondDir);
    clear secondDir
    
    % LOAD STIMULUS PARS AND A DATARUN FOR MASTER
    [StimulusPars,~,~,datarun_mas] = Directories_Params_v23(exp_nm, 'WN', 'mapPRJ');
    SPars.mas = StimulusPars.master;
    SPars.WN  = StimulusPars.slv;
    [StimulusPars] = Directories_Params_v23(exp_nm, 'NSEM', 'mapPRJ');
    SPars.NSEM = StimulusPars.slv;
    clear StimulusPars
    if exist('debug','var') && debug
        SPars.WN.FitBlocks = SPars.WN.FitBlocks(1); 
        SPars.NSEM.FitBlocks = SPars.NSEM.FitBlocks(1); 
    end
    % HACK FIX
    SPars.NSEM.computedtstim = SPars.WN.computedtstim;
    
    
    % LOAD MOVIES
    [blockedmoviecell, inputstats.WN] = loadmoviematfile(exp_nm , 'WN', GLMType.cone_model,'fitmovie');
    concat_fitmovie.WN      = concat_fitmovie_fromblockedcell(blockedmoviecell, SPars.WN);
    clear blockedmoviecell
    [blockedmoviecell, inputstats.NSEM ] = loadmoviematfile(exp_nm , 'NSEM', GLMType.cone_model,'fitmovie');
    concat_fitmovie.NSEM      = concat_fitmovie_fromblockedcell(blockedmoviecell, SPars.NSEM);
    clear blockedmoviecell
    
    
    %%
    for i_celltype = celltypes
        if i_celltype == 1; cellgroup = RCD{i_exp}.ONP;  celltype = 'ONPar'; end
        if i_celltype == 2; cellgroup = RCD{i_exp}.OFFP; celltype = 'OFFPar'; end
        fopt_alltraindata.WN    = cell(length(cellgroup),1);
        fopt_alltraindata.NSEM  = cell(length(cellgroup),1); 
        
       
        
        %%
        for i_cell = 1:length(cellgroup)
            %% 
            % LOAD SPIKE TIMES
            cid = cellgroup(i_cell); 
            cell_savename = sprintf('%s_%d', celltype,cid);
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.WN_organizedspikesdir, cell_savename));
            spikesconcat.WN.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, SPars.WN); 
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.NSEM_organizedspikesdir, cell_savename));
            spikesconcat.NSEM.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, SPars.NSEM);
            
            % LOAD STA AND CENTER COORDINATES 
            eval(sprintf('load %s/STAandROI_%s.mat STAandROI', Dirs.WN_STAdir, cell_savename));
            master_idx         = find(datarun_mas.cell_ids == cid);
            stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
            stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
            slvdim.height      = SPars.WN.height; 
            slvdim.width       = SPars.WN.width; 
            [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, SPars.mas, slvdim);
            cellinfo.WN.WN_STA = STAandROI.STA;
            cellinfo.WN.slave_centercoord = center_coord;
            slvdim.height      = SPars.NSEM.height; 
            slvdim.width       = SPars.NSEM.width; 
            [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, SPars.mas, slvdim);
            cellinfo.NSEM.WN_STA = STAandROI.STA;
            cellinfo.NSEM.slave_centercoord = center_coord;
            clear center_coord sd stafit_centercoord stafit_sd
            
            % LOAD FITTED PARAMETERS FROM RAW_CONV_DATA
            if i_celltype == 1
                P_opt.WN   = RCD{i_exp}.CONV_WN_ONP{i_cell}.fit_p;
                P_opt.NSEM = RCD{i_exp}.CONV_NSEM_ONP{i_cell}.fit_p;
            elseif i_celltype == 2
                P_opt.WN   = RCD{i_exp}.CONV_WN_OFFP{i_cell}.fit_p;
                P_opt.NSEM = RCD{i_exp}.CONV_NSEM_OFFP{i_cell}.fit_p;
            end
            
            %% 
            for i_pct = 1:length(WN_fits)
                [fittedGLM]   = (P_opt.WN spikesconcat.WN, concat_fitmovie.WN, inputstats.WN, cellinfo.WN); 
            end 
            

        end
    end
end


if exist('shortlist','var') && shortlist
    eval(sprintf('save %s/Raw_Conv_shortlist.mat Raw_Conv_Data GLMType', base_savedir))
else
    eval(sprintf('save %s/Raw_Conv.mat Raw_Conv_Data GLMType', base_savedir))
end

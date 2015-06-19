%%
% RASTBPS_WRAP
% AKHEITMAN 2015-03-03: Considered Done
% WRAP of rastbps_comp
% Not STAND ALONE but calls "stand-alone" rastbps_comp
%
% GOAL:
% Compute the log-like of raster rate model with an added post spike term
% Normalization factor for the GLM which includes a conditioning term
% Normalizing with just a pure rate model not always sufficient for GLM
%
%
% INPUT:
% fittedGLM from NSEM_Home/GLM_Output_Analysis directory
% fittedGLM has the recorded raster,fitted ps filters, t_bin
%
% OUTPUT:
% base_crm_findPS_expA.mat file with structure "raster score"
% single structure with raster score for all cells in each experiment
% result put into NSEM_Home/Cell_Selection/Raster_Metrics directory

% AKHEITMAN 2015-03-09 
%    Added plot sections (compute, plotrasterscores, comparetoGLM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRECTORIES AND IDENTIFY WHICH CELLS TO COMPUTE

% RESET
clear; close all; clc

% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE
celltypes = [1 2]
exps = [1] 
debug = false;
compute = true; 
% SET DIRECTORIES / GLM WITH TIME KERNEL
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType);
savedir = sprintf('%s/Raster_Metrics', BD.Cell_Selection);
if ~exist(savedir,'dir'), mkdir(savedir); end

% CRM = CONDITIONALIZED RATE MODEL
%crm_type = 'base_crm_importPS';
%crm_type = 'base_crm_findPS_unitlikebasis';
crm_type = 'base_crm_findPS';
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP LOAD COMPUTE SAVE
if compute
%%% REORGANIZE INTO A MORE APPROACHABLE STRUCTURE
if strcmp(crm_type, 'base_crm_findPS') || strcmp(crm_type, 'base_crm_findPS_HO')
    % HARD PS PARAMETERS  MATCHES GLMPARAMS
    t_bin0 = .00083275;
    ps_params.ms           = 100 ;     
    ps_params.filternumber = 20;
    ps_params.spacing      = pi/2;
    ps_params.bstretch     = .05;
    ps_params.alpha        = 0;
    ps_params.fratio       = .5;  
    ps_basis               = prep_spikefilterbasisGP(ps_params,t_bin0);
end


if strcmp(crm_type, 'base_crm_findPS_unitlikebasis')
    vec_1 = [1;1; zeros(38,1)];
    vec_2 = [1;1;1;1; zeros(76,1)];
    set_1 = NaN(40,20); pad_1 = zeros(80,20);
    set_2 = NaN(80,20); pad_2 = zeros(40,20);
    for i_rotate = 1:20
        rotate_1 = circshift(vec_1,2*(i_rotate-1) );
        rotate_2 = circshift(vec_2,4*(i_rotate-1) );       
        set_1(:,i_rotate) =  rotate_1;
        set_2(:,i_rotate) =  rotate_2;
    end
    
    basis_1 = [set_1; pad_1];
    basis_2 = [pad_2; set_2];
    ps_basis = [basis_1 basis_2];
end


for i_exp = exps
    raster_scores = allcells{i_exp};
    raster_scores.stim_types = {'WN','NSEM'};
    raster_scores.celltypes  = {'ONP','OFFP'};
    
    raster_scores.notes.n1 = 'Both use 10 msec smoothing of raster spike trains';
    raster_scores.notes.n2 = 'UOP: unconditioned optimal performance';
    raster_scores.notes.n3 = 'CRM: conditioned rate model';
    raster_scores.notes.n4 = 'BPS: Bits Per Spike';
    raster_scores.notes.n5 = 'LPPS: logarithmic probability per second';
    raster_scores.notes.n6 = 'GLM: Generalized Linear Model';
    raster_scores.crm_type = crm_type;
    raster_scores.timestamp           = datestr(clock);
    raster_scores.mfile_name          = mfilename('fullpath');
    for i_celltype = celltypes
        % LOAD STIMULUS PARAMETERS / DEFINE CELL NUMBERS
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        map_type= 'mapPRJ';
        if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
        if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
      
        
        if exist('debug','var') && debug
            cellgroup = cellgroup(1:2);
        end
        
        scores_WN.uop_bps    = NaN(length(cellgroup),1);
        scores_WN.crm_bps    = NaN(length(cellgroup),1);
        scores_NSEM.uop_bps  = NaN(length(cellgroup),1);
        scores_NSEM.crm_bps  = NaN(length(cellgroup),1);
        scores_WN.uop_lpps   = NaN(length(cellgroup),1);
        scores_WN.crm_lpps   = NaN(length(cellgroup),1);
        scores_NSEM.uop_lpps = NaN(length(cellgroup),1);
        scores_NSEM.crm_lpps = NaN(length(cellgroup),1);
        
        scores_WN.glm_bps     = NaN(length(cellgroup),1);        
        scores_WN.glm_lpps    = NaN(length(cellgroup),1);
        scores_NSEM.glm_bps   = NaN(length(cellgroup),1);
        scores_NSEM.glm_lpps  = NaN(length(cellgroup),1);
        
        if strcmp(crm_type,'base_crm_findPS')
            scores_WN.psnote1   = 'Raster_Fitted_PSFilter';
            scores_WN.psnote2   = 'Exp of Filter is spike induced gain';
            scores_WN.ps_filter = cell(length(cellgroup),1);
            
            scores_NSEM.psnote1   = 'Raster_Fitted_PSFilter';
            scores_NSEM.psnote2   = 'Exp of Filter is spike induced gain';
            scores_NSEM.ps_filter = cell(length(cellgroup),1);
        end
        cellgroup
        for i_cell = 1:length(cellgroup)
            for i_stimtype = 1:2
                % LOAD FITTED GLM
                if i_stimtype == 1, stimtype = 'WN';   end
                if i_stimtype == 2, stimtype = 'NSEM'; end
                secondDir.exp_nm    = exp_nm;
                secondDir.stim_type = stimtype;
                secondDir.map_type  = 'mapPRJ';
                secondDir.fitname   = GLMType.fitname;
                Dirs.glmfitdir   = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
                clear fittedGLM cell_savename
                cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
                eval(sprintf('load %s/%s.mat', Dirs.glmfitdir, cell_savename))
                
                % PULL OUT COMPONENTS NECESSARY
                inp_rast = fittedGLM.xvalperformance.rasters.recorded;
                t_bin    = fittedGLM.t_bin;
                psfilter = fittedGLM.linearfilters.PostSpike.Filter; 
                
                glm_bps  = fittedGLM.xvalperformance.logprob_glm_bpspike;
                glm_lpps = fittedGLM.xvalperformance.logprob_glm_raw/(t_bin * size(inp_rast,2) ); 
                
                % COMPUTE BITS PER SPIKE FOR CONDITIONED RATE MODEL
                
                   
                if strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    [uop_bits_perspike uop_logprob_persec crm_bits_perspike crm_logprob_persec opt_params] = rastbps_comp_findPS_nosignal(inp_rast,t_bin,ps_basis);                     
                    display(sprintf('%s: Offset=%1.2e, ScaleRate=%1.2e', stimtype,opt_params.offset,opt_params.rate_drive));
                end
                %{
                if strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    [uop_bits_perspike uop_logprob_persec crm_bits_perspike crm_logprob_persec opt_params] = rastbps_comp_findPS(inp_rast,t_bin,ps_basis);                     
                    display(sprintf('%s: Offset=%1.2e, ScaleRate=%1.2e', stimtype,opt_params.offset,opt_params.rate_drive));
                end                
                if strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    [uop_bits_perspike uop_logprob_persec crm_bits_perspike crm_logprob_persec opt_params] =...
                        rastbps_comp_findPS_modsignal(0,inp_rast,t_bin,ps_basis);                     
                    display(sprintf('%s: Offset=%1.2e, ScaleRate=%1.2e', stimtype,opt_params.offset,opt_params.rate_drive));
                end
                %}
                
                
                if i_stimtype == 1 
                    %WN_uop_bps = uop_bits_perspike;
                    %WN_crm_bps = crm_bits_perspike;
                    scores_WN.uop_bps(i_cell)  = uop_bits_perspike;
                    scores_WN.crm_bps(i_cell)  = crm_bits_perspike;
                    scores_WN.uop_lpps(i_cell) = uop_logprob_persec;
                    scores_WN.crm_lpps(i_cell) = crm_logprob_persec;   
                    scores_WN.glm_bps(i_cell)  = glm_bps;
                    scores_WN.glm_lpps(i_cell) = glm_lpps;   
                elseif i_stimtype == 2
                    %NSEM_uop_bps = uop_bits_perspike;
                    %NSEM_crm_bps = crm_bits_perspike;
                    scores_NSEM.uop_bps(i_cell)  = uop_bits_perspike;
                    scores_NSEM.crm_bps(i_cell)  = crm_bits_perspike;
                    scores_NSEM.uop_lpps(i_cell) = uop_logprob_persec;
                    scores_NSEM.crm_lpps(i_cell) = crm_logprob_persec;  
                    scores_NSEM.glm_bps(i_cell)  = glm_bps;
                    scores_NSEM.glm_lpps(i_cell) = glm_lpps;   
                end
                
                if strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    if i_stimtype == 1
                        scores_WN.ps_filter{i_cell}   = opt_params.ps_filter;
                    elseif i_stimtype == 2
                        scores_NSEM.ps_filter{i_cell} = opt_params.ps_filter;
                    end
                end
                
                
            end
            display(sprintf('### %s %s: CRM: LogProb of WN-NSEM: %1.2e ###',...
                expname,cell_savename, (scores_WN.crm_lpps(i_cell)-scores_NSEM.crm_lpps(i_cell))  ));
            %display(sprintf('### %s %s: UOP: LogProb of WN-NSEM: %1.2e ###',...
            %    expname,cell_savename,  (scores_WN.uop_lpps(i_cell)-scores_NSEM.uop_lpps(i_cell)) ));
            
        end
        raster_scores.celltype{i_celltype}.scores_WN = scores_WN;
        raster_scores.celltype{i_celltype}.scores_NSEM = scores_NSEM;
       
        if exist('debug','var') && debug
            eval(sprintf('save %s/%s_%s_nosignal_DBUG.mat raster_scores',savedir,crm_type,exp_nm));
        else
            eval(sprintf('save %s/%s_nosignal_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        end
    end
end

end


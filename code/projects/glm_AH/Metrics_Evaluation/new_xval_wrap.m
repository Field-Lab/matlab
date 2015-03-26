% RECOMPUTE NEW RASTER METRICS
% AKHEITMAN 2014-12-08

%{
clear; clc
i_celltype = 1; i_cell = 1;  i_stimtype = 1;  i_exp = 1;

exps = 1:4; celltypes = 1:2;  stimtypes = 1:2;
metric_type = 'fracvar_avgsignals'
new_xval_wrap(metric_type,exps,stimtypes,celltypes)

clear; clc
exps = 1:4; 
exps = 3
celltypes = 1:2;  stimtypes = 2;
metric_type = 'Viktor_Spike';
new_xval_wrap(metric_type,exps,stimtypes,celltypes)
%}

function new_xval_wrap(metric_type,exps,stimtypes,celltypes)

% JUNK FOR DEBUGGING
%{
clear; clc
exps = 1:4; celltypes = 1:2;  stimtypes = 1:2;
metric_type = 'bits_per_spike'; i_celltype = 1; i_cell = 1;  i_stimtype = 1;  i_exp = 1;
%}
% LOAD CELLS GLOBAL DIRS
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 



% HARD PARAMETERS OF THE RASTER WE ALREADY HAVE SET
hard_params.raster_params.bindur         = .00083275;
hard_params.raster_params.bins_per_frame = 10;
hard_params.map_type = 'mapPRJ';

% SPECIFY GLMTYPE TO GET A FITNAME
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.nullpoint = 'mean';  GLMType.debug = false;
GLMType.CONVEX = true; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.specialchange = false;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.fitname  = GLM_fitname(GLMType); 
fitname = GLMType.fitname



if strcmp(metric_type, 'fracvar_avgsignals')
    rast_met.normalizations = false;
    %rast_met.rawcomputation = true;
    rast_met.foldername     = 'rasterprecision_prediction_avgsignal';
    rast_met.matname_base   = 'rasterprecision_prediction_avgsignal';
    rast_met.params.smoothbins     = [.5 1 2 4 7 10 13 17 20 25 30 40 50 65 80 100 120 150 200 250 350 500 1000];
    rast_met.params.bindur          = .00083275
end

if strcmp(metric_type, 'bits_per_spike')
    rast_met.normalizations         = true;
    rast_met.foldername             = 'rasterprecision_prediction_avgsignal';
    rast_met.matname_base           = 'rasterprecision_prediction_avgsignal';
    rast_met.params.smoothbins      = [.5 1 2 4 7 10 13 17 20 25 30 40 50 65 80 100 120 150 200 250 350 500 1000];
    rast_met.params.bindur          = .00083275
end

if strcmp(metric_type, 'Viktor_Spike')
    rast_met.normalizations         = true;
    rast_met.foldername             = 'rasterprecision_paireddistance_Viktor';
    rast_met.matname_base           = 'rasterprecision_paireddistance_Viktor';
    rast_met.params.ViktorTimeParams_Bins    = [4 8 16 32 64 128 256];
    rast_met.params.bindur                   = hard_params.raster_params.bindur;
end

%%% Values avaliable in the raster metrics %%%
%{

if strcmp(computation_type, 'rasterprecision_prediction_avgsignal')
    metric_params.smoothbins = [.5 1 2 4 7 10 13 17 20 25 30 40 50 65 80 100 120 150 200 250 350 500 1000];
end
if strcmp(computation_type, 'rasterprecision_paireddistance_Viktor')
     metric_params.ViktorTimeParams_Bins = [1 2 4 8 16 32 64 128 256 512 1024];
     metric_params.pairnumbers           = 50;
     metric_params.bindur                = hard_params.raster_params.bindur;    
   % metric_params.ViktorTimeParams_Bins    = [4 8 16 32 64 128 256];% 512 1024];
   % metric_params.pairnumbers	= 30
end
if strcmp(computation_type, 'rasterprecision_paireddistance_Vector')
   metric_params.smoothbins    = [1 2 4 8 16 32 64 128 256 512 1024];
   metric_params.pairnumbers	= 50; 
end
%}
%%
for i_exp = exps
    for i_stimtype = stimtypes  
        for i_celltype = celltypes
            %% Experiment Dependent Parameters
            

            % LOAD CELLS
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
           
            % RASTER METRIC DIRECTORY (FOR NORMALIZATION)
            secondDir.exp_nm     = exp_nm;
            secondDir.stim_type  = stimtype;
            secondDir.map_type   = hard_params.map_type;
            Dirs.organizedspikes = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
            Dirs.rastmet         = sprintf('%s/%s', Dirs.organizedspikes, rast_met.foldername);
            
            % FITTED GLM
            secondDir.fitname   = fitname;
            Dirs.fittedGLM      = NSEM_secondaryDirectories('savedir_GLMfit', secondDir);  clear secondDir
            
            % SAVE DIR 
            Dirs.crossval_savedir = sprintf('%s/crossval_%s', Dirs.fittedGLM, metric_type);
            if ~exist(Dirs.crossval_savedir, 'dir'), mkdir(Dirs.crossval_savedir); end
            
            cellgroup = fliplr(cellgroup)
            
            %%
            for i_cell = 1:length(cellgroup)
                % CLEAN UP 
                clear crossval_perf scores
                
                % LOAD STORED RASTER METRICS
                cid = cellgroup(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
                display(sprintf('Working on Exp: %s Stim: %s Cell: %s', exp_nm,stimtype, cell_savename));
                
                % LOAD FITTEDGLM ->EXTRACT RASTERS
                eval(sprintf('load %s/%s.mat', Dirs.fittedGLM, cell_savename));
                
                rast_sim = fittedGLM.xvalperformance.rasters.glm_sim;
                rast_rec = fittedGLM.xvalperformance.rasters.recorded;
                
                
                
                crossval_perf.cid                  = cid;
                crossval_perf.exp_nm               = exp_nm;
                crossval_perf.celltype             = celltype;
                crossval_perf.name                 = cell_savename;
                crossval_perf.fitname              = fitname;
                crossval_perf.metric_type          = metric_type;
                crossval_perf.metric_pars          = rast_met;
                            
                
                % COMPUTE RAW SCORES
                if strcmp(metric_type,'fracvar_avgsignals')
                    [scores] = rawcrossval_performance_fracvar(rast_sim,rast_rec,rast_met.params.smoothbins);
                elseif strcmp(metric_type, 'bits_per_spike')
                    scores.metric_raw          = fittedGLM.xvalperformance.logprob_glm_bpspike *ones(size(rast_met.params.smoothbins));  
                    scores.notes.bps_eqn       = 'LOGPROB(RAST|MODEL) - LOGPROB(RAST|MEAN RATE)';
                    scores.notes.bps_time      = 'RAW SCORE IS TIME INDEPENDENT';
                    scores.notes.normfactor    = 'BPS of the raster avg signal (avg signal determined at different time scales';
                    scores.sigma_bin_norm      = rast_met.params.smoothbins;
                elseif strcmp(metric_type, 'Viktor_Spike')
                    [scores] = rawcrossval_performance_ViktorSpike(rast_sim,rast_rec,rast_met.params.ViktorTimeParams_Bins,rast_met.params.bindur);
                end
                
                % NORMALIZE IF NORMALIZABLE 
                if rast_met.normalizations
                    
                    % LOAD RASTERMETRICS
                    matfilename   = sprintf('%s_%s', rast_met.matname_base, cell_savename);
                    eval(sprintf('load %s/%s.mat', Dirs.rastmet, matfilename));
                    
                    % 
                    if strcmp(metric_type, 'bits_per_spike')
                        normfactor = raster_metrics.prediction_avgsignal.bps_mean';
                    elseif strcmp(metric_type, 'Viktor_Spike')
                        rast_viktor =  raster_metrics.paireddistances_viktor.Viktordist_perspike;
                        rast_bins   =  raster_metrics.paireddistances_viktor.ViktorParams_Bins;
                        [dummy , ind_rast ] = intersect(rast_bins,rast_met.params.ViktorTimeParams_Bins);
                        
                        normfactor = rast_viktor(ind_rast);
                        clear dummy ind_rast rast_bins rast_viktor
                    end
                    scores.metric_normsubtract = scores.metric_raw - normfactor;
                    scores.metric_normdivide   = scores.metric_raw./normfactor;
                    
                end
                crossval_perf.scores = scores;
                
                
                % SAVE
                eval(sprintf('save %s/crossvalperf_%s.mat crossval_perf', Dirs.crossval_savedir, cell_savename))
            end


        
        end
    end
end


end
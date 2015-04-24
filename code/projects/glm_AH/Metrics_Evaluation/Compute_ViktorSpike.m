% AKHEITMAN 2015-01-19
% Minimal Version with 13 and 32 Bins for ViktorTimeParams_Bins


%{
clear
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
changes_cell{2}.type = 'input_pt_nonlinearity';
changes_cell{2}.name = 'piecelinear_fourpiece_eightlevels';
exps = [1 2 3 4];
Compute_ViktorSpike(exps,changes_cell) ;


clear
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
changes_cell{2}.type = 'input_pt_nonlinearity';
changes_cell{2}.name = 'piece_linear_aboutmean';
exps = [1 2 3 4];
Compute_ViktorSpike(exps,changes_cell) ;


clear
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_fullcone';
exps = [1 2 3 4];
Compute_ViktorSpike(exps,changes_cell) ;

clear
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'rk2-ConductanceBased';
exps = [3];
Compute_ViktorSpike(exps,changes_cell) ;


exps = [1 2 3 4];
Compute_ViktorSpike(exps) ;

%}
%%
function Compute_ViktorSpike(exps,changes_cell) 

BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 


stimtypes =  [1 2];
celltypes = [1 2];

% HARD PARAMETERS OF THE RASTER WE ALREADY HAVE SET
hard_params.raster_params.bindur         = .00083275;
hard_params.raster_params.bins_per_frame = 10;
hard_params.map_type = 'mapPRJ';

metric_type = 'Viktor_Spike';

%changes_cell{1}.type = 'filter_mode';
%changes_cell{1}.name = 'fixedSP-ConductanceBased';
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType); 
GLMType.func_sname = 'glmwrap';
GLMType.fullmfilename =mfilename('fullpath'); 
fitname = GLMType.fitname;

cellselectiontype = 'shortlist';
if strcmp(metric_type, 'Viktor_Spike')
    rast_met.normalizations         = true;
    rast_met.foldername             = 'rasterprecision_paireddistance_Viktor';
    rast_met.matname_base           = 'rasterprecision_paireddistance_Viktor';
    rast_met.params.ViktorTimeParams_Bins    = [13 32];  % 10 msec and 25 msec 
    rast_met.params.bindur                   = hard_params.raster_params.bindur;
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
           
            if strcmp(cellselectiontype, 'shortlist')
                [exp_nm,cells,expname]  = cell_list( i_exp, 'shortlist'); cells = cell2mat(cells);
                cellgroup = intersect(cellgroup, cells);
            end
            
            % RASTER METRIC DIRECTORY (FOR NORMALIZATION)
            secondDir.exp_nm     = exp_nm;
            secondDir.stim_type  = stimtype;
            secondDir.map_type   = 'mapPRJ';

            % FITTED GLM
            secondDir.fitname   = fitname;
            Dirs.fittedGLM      = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);  clear secondDir
            
            % SAVE DIR 
            Dirs.crossval_savedir = sprintf('%s/crossval_%s', Dirs.fittedGLM, metric_type);
            if ~exist(Dirs.crossval_savedir, 'dir'), mkdir(Dirs.crossval_savedir); end            
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
                
                crossval_perf.scores = scores;
                
                
                % SAVE
                eval(sprintf('save %s/crossvalperf_%s.mat crossval_perf', Dirs.crossval_savedir, cell_savename))
            end


        
        end
    end
end


end

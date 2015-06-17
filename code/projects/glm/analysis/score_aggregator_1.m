% GOAL:
% Have a way to continually update and collect scores
% Either Creates or Updates the Aggregated Scores structure
% Outputs to Analysis Directory before the WN and NSEM fork
% Version 0 works.. not sure how auto-update is going though 2015-06-16
% Version 1 has cleaner auto-update and integrated viktor spike

% AKHEITMAN 2015-06-16
% CALLS:
% GLM_settings, GLM_fitname NSEM_BaseDirectories cell_list


% Standard Bookkeeping
clear ; close all; clc;

% Dictate Model
%model.settings = {};
model.settings{1}.type = 'PostSpikeFilter';
model.settings{1}.name =  'OFF';
%model.settings{1}.type = 'PostSpikeFilter';
%model.settings{1}.name =  'OFF';
%model.settings{1}.type = 'cone_model';
%model.settings{1}.name = 'rieke_linear';
%model.settings{1}.type = 'cone_model';
%model.settings{1}.name = 'rieke_fullcone';

%details.metric  = 'crossval_BPS';
%details.metric_note = 'Bits Per Spike over crossvalidated dataset: (logprob(rast|model)-logprob(rast|flatrate))/spikes';
%details.metric      = 'crossval_fracvar_10msec';
%details.metric_note  = 'Fraction of Variance Explained: CrossValidated Dataset';

%details.metric       = 'crossval_victorspike_50msec';
%details.metric_note  = 'Victor Spike with 50 msec timescale: CrossValidated Dataset';


details.metric       = 'crossval_victorspike_10msec'
details.metric_note  = 'Victor Spike with 10 msec timescale: CrossValidated Dataset';

details.exps     = [1 2 3 4];

% LOAD DATA FROM FITTED GLM
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
model.GLMType     = GLM_settings('default',model.settings);
model.fitname     = GLM_fitname(model.GLMType);
savedir           = sprintf('%s/%s', BD.GLM_output_analysis, model.fitname);




for i_exp = details.exps
    exp_nm  = allcells{i_exp}.exp_nm;
    filename = sprintf('%s/%s_%s.mat', savedir,details.metric, exp_nm);
    expname = allcells{i_exp}.expname;
    
    
    % Create the aggregated scores structure
    if ~exist(filename)
        aggregated_scores = allcells{i_exp};  
        aggregated_scores.metric = details.metric;
        aggregated_scores.metric_note = details.metric_note;
        for i_celltype = [1,2]
            if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
            if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
            aggregated_scores.celltype{i_celltype}.scores_WN   = NaN(length(cellgroup),1);
            aggregated_scores.celltype{i_celltype}.scores_NSEM = NaN(length(cellgroup),1);
        end   
    else
        display(sprintf('checking for updates to %s on %s', model.fitname,exp_nm));
        eval(sprintf('load %s/%s_%s.mat aggregated_scores', savedir,details.metric, exp_nm));
        combined_scores = [aggregated_scores.celltype{1}.scores_WN; aggregated_scores.celltype{1}.scores_NSEM;...
           aggregated_scores.celltype{2}.scores_WN; aggregated_scores.celltype{1}.scores_NSEM];
        exist_nan = sum(isnan(combined_scores));
        if exist_nan == 0
           display('No Updates Necessary .. All scores Present')
           continue
        end       
    end
    
    for i_celltype = [1,2]
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        
        for i_fit = 1:2
            % Pull the scores
            if i_fit == 1, rawscores = aggregated_scores.celltype{i_celltype}.scores_WN; end
            if i_fit == 2, rawscores = aggregated_scores.celltype{i_celltype}.scores_NSEM; end
            
            % Find the correct load directory
            secondDir.exp_nm    = exp_nm;
            secondDir.fitname = model.fitname;
            secondDir.map_type  = 'mapPRJ';
            if i_fit == 1, secondDir.stim_type = 'WN'; end
            if i_fit == 2, secondDir.stim_type = 'NSEM'; end  
            loaddir = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
            clear secondDir

            stillNAN = find(isnan(rawscores));    
            % Run through cells that still have NaN as the entry
            for i_cell = 1:length(stillNAN);
                i_index = stillNAN(i_cell);
                clear crossval_rawmetrics cell_savename
                cid = cellgroup(i_index);
                cell_savename = sprintf('%s_%d', celltype,cid);
                matfilename = sprintf('%s/%s.mat', loaddir, cell_savename);
                
                % Compute metric if necessary
                if exist(matfilename)
                    display(sprintf('LOADING %s %s', expname,cell_savename));
                    eval(sprintf('load %s fittedGLM', matfilename));

                    if strcmp(details.metric,'crossval_BPS')
                        rawscores(i_index) = fittedGLM.xvalperformance.logprob_glm_bpspike;
                    end
                    if strcmp(details.metric,'crossval_fracvar_10msec')
                        smoothbins = 12;
                        rast_rec = fittedGLM.xvalperformance.rasters.recorded;
                        rast_sim = fittedGLM.xvalperformance.rasters.glm_sim; 
                        [crossval_score] = rawcrossval_performance_fracvar(rast_sim,rast_rec,smoothbins);
                        rawscores(i_index) = crossval_score.metric_raw;
                    end
                    if strcmp(details.metric,'crossval_victorspike_50msec')
                        timescale = .05;
                        bindur = fittedGLM.t_bin;
                        rast_rec = fittedGLM.xvalperformance.rasters.recorded;
                        rast_sim = fittedGLM.xvalperformance.rasters.glm_sim; 
                        [crossval_score] = rawcrossval_performance_ViktorSpike(rast_sim,rast_rec,timescale, bindur);
                        rawscores(i_index) = crossval_score.metric_raw;
                    end
                    if strcmp(details.metric,'crossval_victorspike_10msec')
                        timescale = .01;
                        bindur = fittedGLM.t_bin;
                        rast_rec = fittedGLM.xvalperformance.rasters.recorded;
                        rast_sim = fittedGLM.xvalperformance.rasters.glm_sim; 
                        [crossval_score] = rawcrossval_performance_ViktorSpike(rast_sim,rast_rec,timescale, bindur);
                        rawscores(i_index) = crossval_score.metric_raw;
                    end
                    
                end
            end

            if i_fit == 1
                aggregated_scores.celltype{i_celltype}.scores_WN = rawscores;
            end
            if i_fit == 2
                aggregated_scores.celltype{i_celltype}.scores_NSEM = rawscores;
            end                
        end
    end
    aggregated_scores.timestamp       = datestr(clock);
    aggregated_scores.code_name       = mfilename('fullpath');
    eval(sprintf('save %s/%s_%s.mat aggregated_scores', savedir,details.metric, exp_nm))
end 



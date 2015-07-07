%{
clear ; close all; clc;
special_arg = '';
metric_type.name  = 'fracvar_10msec_ODDEVEN';
metric_type.note  = 'Fraction of Variance Explained: Rates from ODD and EVEN Portions of Raster';
exps     = [1 2 3 4];
normalizer_aggregator(metric_type,exps,special_arg)
%}
function normalizer_aggregator(metric_type,exps, special_arg)
% Standard Bookkeeping
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
savedir = sprintf('%s/Raster_Metrics', BD.Cell_Selection);


% Hacked way of loading the recorded raster (no sim!)
glm_settings = {};
GLMType     = GLM_settings('default',glm_settings);
fitname     = GLM_fitname(GLMType);
clear glm_settings GLMType





for i_exp = exps
    exp_nm  = allcells{i_exp}.exp_nm;
    filename = sprintf('%s/%s_%s.mat', savedir,metric_type.name, exp_nm)
    expname = allcells{i_exp}.expname;
    
    raster_scores = allcells{i_exp};  
	raster_scores.metric = metric_type.name;
	raster_scores.metric_note = metric_type.note;
	for i_celltype = [1,2]
            if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
            if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
            raster_scores.celltype{i_celltype}.scores_WN   = NaN(length(cellgroup),1);
            raster_scores.celltype{i_celltype}.scores_NSEM = NaN(length(cellgroup),1);
	end
    % Create the aggregated scores structure
    for i_celltype = [1,2]
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        
        for i_fit = 1:2
            % Pull the scores
            if i_fit == 1, rawscores = raster_scores.celltype{i_celltype}.scores_WN; end
            if i_fit == 2, rawscores = raster_scores.celltype{i_celltype}.scores_NSEM; end
            
            % Find the correct load directory
            secondDir.exp_nm    = exp_nm;
            secondDir.fitname   = fitname;
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
                    eval(sprintf('load %s', matfilename));

                    if strcmp(metric_type.name,'fracvar_10msec_ODDEVEN')
                        smoothbins = 12;
                        rast_rec = fittedGLM.xvalperformance.rasters.recorded;
                        
                        reps = size(rast_rec,1);
                        if mod(reps,2) == 1
                            reps = reps-1;
                        end
                        odd_reps = [1:2:reps];
                        even_reps = [2:2:reps];
                        
                        rast_odd = rast_rec(odd_reps,:);
                        rast_even = rast_rec(even_reps,:);

                        [crossval_score] = rawcrossval_performance_fracvar(rast_odd,rast_even,smoothbins);
                        rawscores(i_index) = crossval_score.metric_raw;
                    end
                end
            end

            if i_fit == 1
                raster_scores.celltype{i_celltype}.scores_WN = rawscores;
            end
            if i_fit == 2
                raster_scores.celltype{i_celltype}.scores_NSEM = rawscores;
            end
        end
    end
    raster_scores.timestamp       = datestr(clock);
    raster_scores.code_name       = mfilename('fullpath');
    eval(sprintf('save %s/%s_%s.mat raster_scores', savedir,metric_type.name, exp_nm))
end


end




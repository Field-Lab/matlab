for i_exp = exptests
    % LOAD CELLS  (OLD WAY)
    expnumber = i_exp;
    [exp_nm] = allcells{i_exp}.exp_nm;
    
    
    %%%%%% Identify Directories %%%%%%%%%%%
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = 'mapPRJ'; 
    inputs.stim_type = 'WN';
    inputs.fitname   = GLMType.fitname;
    d_load_WN = NSEM_secondaryDirectories('loaddir_GLMfit', inputs);  clear inputs; 
    
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = 'mapPRJ'; 
    inputs.stim_type = 'NSEM';
    inputs.fitname   = GLMType.fitname;
    d_load_NSEM = NSEM_secondaryDirectories('loaddir_GLMfit', inputs);  clear inputs; 
    

    for i_metric = 1:4
        model_comparison.scores{i_exp}.metricnames{i_metric} =  metrics.plots.title_base{i_metric};
    end
    model_comparison.scores{i_exp}.note = 'Column 1 is WN, Column 2 is NSEM';
%%
    for i_celltype = celltypes
        
        if i_celltype == 1
            cell_list = allcells{i_exp}.ONP;
        elseif i_celltype == 2
            cell_list = allcells{i_exp}.OFFP;
        end
%%
        for i_cell = 1:length(cell_list)
            cid = cell_list(i_cell);
            if i_celltype == 1, cell_savename  = sprintf('ONPar_%d',cid);       end
            if i_celltype == 2, cell_savename  = sprintf('OFFPar_%d',cid);      end
            if i_celltype == 1, rast_normindex = find(raster_scores.ONP == cid); end
            if i_celltype == 2, rast_normindex = find(raster_scores.OFFP== cid); end
            
            eval(sprintf('load %s/%s.mat', d_load_WN,cell_savename));
            xval_WN = fittedGLM.xvalperformance; clear fittedGLM
            eval(sprintf('load %s/%s.mat', d_load_NSEM,cell_savename));
            xval_NSEM = fittedGLM.xvalperformance; clear fittedGLM


            score_WN_raw    =    xval_WN.logprob_glm_bpspike; 
            score_NSEM_raw  =  xval_NSEM.logprob_glm_bpspike;
            normalizer_WN_uop   = raster_scores.stim_type{1}.celltype{i_celltype}.scores.values(rast_normindex,1);
            normalizer_NSEM_uop = raster_scores.stim_type{2}.celltype{i_celltype}.scores.values(rast_normindex,1);
            normalizer_WN_cop   = raster_scores.stim_type{1}.celltype{i_celltype}.scores.values(rast_normindex,2);
            normalizer_NSEM_cop = raster_scores.stim_type{2}.celltype{i_celltype}.scores.values(rast_normindex,2);
  
          
        end
        
    end
    
end

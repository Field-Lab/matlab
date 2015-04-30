% AKHEITMAN 2015-03-01
% TRY TO COME UP WITH "CORRECT NORMALIZATION FACTOR"
% COP TESTING  (CONDINTIONALIZED OPTIMUM)

clear; close all; clear all; clc


BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire

exptests = [1 2 3 4];
celltypes = [1 2];
baseoutput_dir = '/Users/akheitman/NSEM_Home/PrototypePlots/Performance_Comparisons/WNWN_vs_NSEMNSEM';
datainput_dir  = '/Users/akheitman/NSEM_Home/PrototypePlots/input_data';
if ~exist(baseoutput_dir), mkdir(baseoutput_dir); end
cellselectiontype = 'shortlist';
cellselectiontype = 'all';
%celltypes = 2; celltype_string = 'OFF';

metrics.plots.init_minval = Inf;
metrics.plots.init_maxval = -Inf;
metrics.plots.title_base{1}  = 'BPS(GLM) - BPS(Optimal Rate Model(10msec))' ;
metrics.plots.title_base{2}  = 'BPS(GLM) - BPS(Linear Conditioned Model(10msec))' ;
metrics.plots.title_base{3}  = 'BPS(GLM) / BPS(Optimal Rate Model(10msec))' ;
metrics.plots.title_base{4}  = 'BPS(GLM) / BPS(Linear Conditioned Model(10msec))' ;
metrics.plots.xlabel         = 'White Noise';
metrics.plots.ylabel         = 'Natural Scenes';


% OLD WAY OF LOADING %
% eval(sprintf('load %s/%s',datainput_dir,metrics.raster_normalization_structure))
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType); 
GLMType.func_sname = 'glmwrap';
GLMType.fullmfilename =mfilename('fullpath'); 
%metrics.shortname = 'fournormPSHO';
metrics.shortname = 'base_crm_importPS';

expstring = '';
for i_exp = exptests;
    expstring = sprintf('%s%d',expstring, i_exp);
end




model_comparison.metrics        = metrics;
model_comparison.fitname        = GLMType.fitname;
model_comparison.fullGLMType    = GLMType;
model_comparison.note1          = 'rows are experiments, columns are the different models';
model_comparison.note2          = 'row1 is WN scores, row2 is NSEM scores';



%%
for i_exp = exptests
    % LOAD CELLS  (OLD WAY)
    expnumber = i_exp
    [exp_nm] = allcells{i_exp}.exp_nm;
    
    %{
    if i_exp == 1,load /Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics/base_crm_findPS_expA.mat; end
    if i_exp == 2,load /Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics/base_crm_findPS_expB.mat; end
    if i_exp == 3,load /Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics/base_crm_findPS_expC.mat; end
    if i_exp == 4,load /Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics/base_crm_expD.mat; end
    %}
    metdir = '/Users/akheitman/NSEM_Home/Cell_Selection/Raster_Metrics';
    eval(sprintf('load %s/base_crm_importPS_%s', metdir, exp_nm));
    
   % if strcmp(cellselectiontype, 'shortlist')
   %     [~,cell_subset,expname]  = cell_list( expnumber, 'shortlist'); cells = cell2mat(cells);        
   % end
    
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
            cells = allcells{i_exp}.ONP;
            if exist('cell_subset','var')
                cells = intersect(cell_subset, cells);
            end
            model_comparison.scores{i_exp}.ONP_cells  = cells;
            
            for i_metric = 1:4                
                model_comparison.scores{i_exp}.ONP_scores{i_metric} = NaN(2,length(cells));
            end
        elseif i_celltype == 2
            cells = allcells{i_exp}.OFFP;
            if exist('cell_subset','var')
                cells = intersect(cell_subset, cells);
            end
            model_comparison.scores{i_exp}.OFFP_cells  = cells;
            for i_metric = 1:4
                model_comparison.scores{i_exp}.OFFP_scores{i_metric} = NaN(2,length(cells));
            end
        end
%%
        for i_cell = 1:length(cells)
            cid = cells(i_cell);
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
            
            for i_metric = 1:4
                if i_metric == 1, score_WN = score_WN_raw - normalizer_WN_uop; end
                if i_metric == 2, score_WN = score_WN_raw - normalizer_WN_cop; end
                if i_metric == 3, score_WN = score_WN_raw / normalizer_WN_uop; end
                if i_metric == 4, score_WN = score_WN_raw / normalizer_WN_cop; end
                
                if i_metric == 1, score_NSEM = score_NSEM_raw - normalizer_NSEM_uop; end
                if i_metric == 2, score_NSEM = score_NSEM_raw - normalizer_NSEM_cop; end
                if i_metric == 3, score_NSEM = score_NSEM_raw / normalizer_NSEM_uop; end
                if i_metric == 4, score_NSEM = score_NSEM_raw / normalizer_NSEM_cop; end
                
                if i_celltype == 1
                    model_comparison.scores{i_exp}.ONP_scores{i_metric}(1,i_cell)  = score_WN;
                    model_comparison.scores{i_exp}.ONP_scores{i_metric}(2,i_cell)  = score_NSEM; 
                elseif i_celltype == 2
                    model_comparison.scores{i_exp}.OFFP_scores{i_metric}(1,i_cell) = score_WN;
                    model_comparison.scores{i_exp}.OFFP_scores{i_metric}(2,i_cell) = score_NSEM;                   
                end
            end  
          
        end
        
    end
    
end




%% Plotting WNWN and NSEMNSEM New
savedir  = sprintf ('%s/%s/',  baseoutput_dir, GLMType.fitname);
if ~exist(savedir, 'dir')
    mkdir(savedir);
end
% Square Plots
minval = metrics.plots.init_minval;
maxval = metrics.plots.init_maxval;

MS_A = 16;
MS_B = 20;


BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection));  % allcells structire

conv_criterion = allcells_glmconv{1}.conv_criterion;

%%
for i_conv = 0:4
	
    for i_label = 1
        if i_label == 1, label_points = false; end
        if i_label == 2, label_points = true; end
    end
    if i_conv == 0
        convcut = 0;
    else
        convcut = conv_criterion(i_conv);
    end
    %figure(1); clf; hold on;
    for i_exp = exptests
        if i_exp == 1, colorstring = 'r'; end
        if i_exp == 2, colorstring = 'g'; end
        if i_exp == 3, colorstring = 'b'; end
        if i_exp == 4, colorstring = 'c'; end
        
        for i_celltype = celltypes

            
            if i_celltype == 1
                celltype_string = 'ONP'; 
                scores = model_comparison.scores{i_exp}.ONP_scores; cids = model_comparison.scores{i_exp}.ONP_cells;
                
                if i_conv == 0
                    passindex = ones(length(cids),1);
                else
                    passindex = find(allcells_glmconv{i_exp}.ONP_CONV(:,i_conv));
                end
                %scores = scores_raw(passindex);
            end
            
             if i_celltype == 2
                celltype_string = 'OFFP'; 
                scores = model_comparison.scores{i_exp}.OFFP_scores; cids = model_comparison.scores{i_exp}.OFFP_cells;
                
                if i_conv == 0
                    passindex = ones(length(cids),1);
                else
                    passindex = find(allcells_glmconv{i_exp}.OFFP_CONV(:,i_conv));
                end
                %scores = scores_raw(passindex);
            end
                

            savename = sprintf('%s_%1.1e_%s_%s_%s',  metrics.shortname,convcut,expstring, celltype_string, cellselectiontype);
            
            for i_metric = 1:4
                subplot(2,2,i_metric); hold on; set(gca,'fontsize', 10)
                vals = scores{i_metric};
                vals = vals(:,passindex);
                
                minval = Inf;
                maxval = -Inf;
                bad_WN   = union(find(abs(vals(1,:))==Inf), find(vals(1,:)==NaN));
                vals(:,bad_WN) = [];
                cids(:,bad_WN) = [];

                bad_NSEM = union(find(abs(vals(2,:)) == Inf), find(vals(2,:)== NaN));
                vals(:,bad_NSEM) = [];
                cids(:,bad_NSEM) = [];


                if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
                if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
                plotstring_center = sprintf('%s.',colorstring);

                if strcmp(transform,'ID')
                    vals = vals;
                    vals = vals;
                end
                if length(transform) >=8 &&  strcmp(transform(1:8),'logistic')
                    vals = 1./ ( 1 + exp(-(log_k)*vals) ) - .5;
                    vals = 1./ ( 1 + exp(-(log_k)*vals) ) - .5;
                end
                
                newmax = max(vals(:));
                dummy  = sort(vals(:));
                newmin = dummy(5);
                if i_metric ==3 || i_metric ==4
                    newmin = 0;
                end
                
                vals(find(vals<=newmin)) = newmin;
                maxval = max( maxval, newmax );
                minval = min( minval, newmin );
                
                if i_metric ==3 || i_metric == 4
                    maxval = max(1,maxval)
                end
                % for all cells knock out bottom 5 -10 values
                plot(vals(1,:), vals(2,:), plotstring_center, 'markersize', MS_B);
                plot(vals(1,:), vals(2,:), plotstring_surround, 'markersize', MS_A);

                if exist('label_points','var') && label_points
                    for i_cell = 1:length(cids)
                        lab = sprintf('%d', cids(i_cell));
                        text(vals(1,i_cell),vals(2,i_cell),lab);
                    end
                end
                unity_line = linspace(minval, maxval,100);
                plot(unity_line,unity_line,'k')
                xlim([minval,maxval]);
                ylim([minval,maxval]);
                xlabel(metrics.plots.xlabel);
                ylabel(metrics.plots.ylabel);
                
                
                title(sprintf('Metric: %s',metrics.plots.title_base{i_metric}),'interpreter','none') 
                
                plotname = sprintf('%s', transform);

                if exist('label_points','var') && label_points
                    plotname = sprintf('%s_Labels',plotname);
                end
                orient landscape
                eval(sprintf('print -dpdf %s/%s_%s.pdf', savedir, savename, plotname));
            end
                
                
                
                %{
                figure(1); subplot(2,2,i_metric); hold on; set(gca,'fontsize', 10)
                % for all cells knock out bottom 5 -10 values
                plot(vals(1,:), vals(2,:), plotstring_center, 'markersize', MS_B);
                plot(vals(1,:), vals(2,:), plotstring_surround, 'markersize', MS_A);

                if exist('label_points','var') && label_points
                    for i_cell = 1:length(cids)
                        lab = sprintf('%d', cids(i_cell));
                        text(vals(1,i_cell),vals(2,i_cell),lab);
                    end
                end
                unity_line = linspace(minval, maxval,100);
                plot(unity_line,unity_line,'k')
                if i_metric == 3 || i_metric ==4
                    xlim([minval,maxval]);
                    ylim([minval,maxval]);
                end
                xlabel(metrics.plots.xlabel);
                ylabel(metrics.plots.ylabel);
                
                
                title(sprintf('Metric: %s',metrics.plots.title_base{i_metric}),'interpreter','none') 
                
                plotname = sprintf('%s', transform);

                if exist('label_points','var') && label_points
                    plotname = sprintf('%s_Labels',plotname);
                end
                %}
            
        end
    end
    %{
    figure(1)
    orient landscape
    eval(sprintf('print -dpdf %s/%1.1e_ALL_%s.pdf', savedir,convcut, savename, plotname));
    %}

end


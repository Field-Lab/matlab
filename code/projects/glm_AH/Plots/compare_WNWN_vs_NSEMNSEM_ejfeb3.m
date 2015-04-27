% Specifically tuned to give the plots that EJ wants for 2015, Feb 3
% started 1-31-2015

% AKHEITMAN 2015-01-26
% USE THE DECIDED METRICS TO MAKE DEFINITIVE WHITE NOISE BEATS NSEM PLOTS
% VKSP with 25 msec (32 bins) as the time scale of cost parameter
% BitsPerSpike with Normalized by rate modeled smoothed with sigma = 10 msec

% Output: NSEM_Home/PerformanceComparisons
% Works on 2015-01-25

clear; close all; clear all; clc
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire

baseoutput_dir = '/Users/akheitman/NSEM_Home/PrototypePlots/Performance_Comparisons/WNWN_vs_NSEMNSEM';
datainput_dir  = '/Users/akheitman/NSEM_Home/PrototypePlots/input_data';
if ~exist(baseoutput_dir), mkdir(baseoutput_dir); end
cellselectiontype = 'shortlist';
exptests = [1 2 3];

%changes_cell{1}.type = 'cone_model';
%changes_cell{1}.name = 'rieke_fullcone';  shortname = 'Dynamic Cone Model';
%changes_cell{1}.type = 'filter_mode';
%changes_cell{1}.name = 'fixedSP-ConductanceBased';
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear'; %shortname = 'Standard GLM';
changes_cell{2}.type = 'input_pt_nonlinearity';
changes_cell{2}.name = 'piecelinear_fourpiece_eightlevels'; shortname = '"General" Input-NL: 4 Piece-Linear';
%changes_cell{2}.type = 'input_pt_nonlinearity';
%changes_cell{2}.name = 'piece_linear_aboutmean'; shortname = 'Hinge About Mean: One Param NL';

% LOOP THROUGH BOTH METRICS
for i_metric = 1
    clear metrics
    if i_metric == 1, metrics.name ='BPS_normsubtract_smooth10msec'; metrics.shortname = 'BPS'; end
    if i_metric == 2, metrics.name ='Viktor_normsubtract_25msec'; metrics.shortname = 'VKSP'; end

if strcmp(metrics.name,'BPS_normsubtract_smooth10msec')
    metrics.raster_normalization_structure = 'selfprediction_10msecsmooth_bps.mat';
    metrics.plots.init_minval = 0;
    metrics.plots.init_maxval = 0;
    metrics.plots.title_base  = 'BPS(GLM)-BPS(Optimal Rate Model(10msec))'
    metrics.plots.xlabel      = ' White Noise';
    metrics.plots.ylabel      = 'Natural Scenes';
end
if strcmp(metrics.name,'Viktor_normsubtract_25msec')
    metrics.raster_normalization_structure = 'rasterprecision_paireddistance_Viktor_25msec.mat';
    metrics.plots.init_minval = 0;
    metrics.plots.init_maxval = 1;
    metrics.plots.title_base  = 'VSP(Sim to Raster) - VSP(Raster to Raster)'
    metrics.plots.xlabel      = ' White Noise';
    metrics.plots.ylabel      = 'Natural Scenes';
end



eval(sprintf('load %s/%s',datainput_dir,metrics.raster_normalization_structure))



if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType); 
GLMType.func_sname = 'glmwrap';
GLMType.fullmfilename =mfilename('fullpath'); 


savedir  = sprintf ('%s/%s/',  baseoutput_dir, GLMType.fitname);
savename = sprintf('%s_%s',  metrics.shortname,cellselectiontype);

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

model_comparison.metrics        = metrics;
model_comparison.fitname        = GLMType.fitname;
model_comparison.fullGLMType    = GLMType;
model_comparison.note1          = 'rows are experiments, columns are the different models';
model_comparison.note2          = 'row1 is WN scores, row2 is NSEM scores';

i_exp = 1;
i_celltype = 1;
i_cell = 1;


%%
for i_exp = exptests
    
    % LOAD CELLS  (OLD WAY)
    expnumber = i_exp;
    if strcmp(cellselectiontype, 'shortlist')
        [exp_nm,cells,expname]  = cell_list( expnumber, 'shortlist'); cells = cell2mat(cells);
    end
    
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
    
    if strcmp(metrics.name,'Viktor_normsubtract_25msec')
        d_load_WN   = sprintf('%s/crossval_Viktor_Spike', d_load_WN);
        d_load_NSEM = sprintf('%s/crossval_Viktor_Spike', d_load_NSEM);
    end
%%
    for i_celltype = 1:2
        
        if i_celltype == 1
            celllist = intersect(cells, allcells{i_exp}.ONP);
            model_comparison.scores{i_exp}.ONP_cells  = celllist;
            model_comparison.scores{i_exp}.ONP_scores = NaN(2,length(celllist));
        elseif i_celltype == 2
            celllist = intersect(cells, allcells{i_exp}.OFFP);
            model_comparison.scores{i_exp}.OFFP_cells  = celllist;
            model_comparison.scores{i_exp}.OFFP_scores = NaN(2,length(celllist));
        end

        for i_cell = 1:length(celllist)
            cid = celllist(i_cell);
            if i_celltype == 1, cell_savename  = sprintf('ONPar_%d',cid);       end
            if i_celltype == 2, cell_savename  = sprintf('OFFPar_%d',cid);      end
            if i_celltype == 1, rast_normindex = find(raster_scores{i_exp}.ONP == cid); end
            if i_celltype == 2, rast_normindex = find(raster_scores{i_exp}.OFFP== cid); end
            
            if strcmp(metrics.name,'Viktor_normsubtract_25msec')
                cell_savename = sprintf('crossvalperf_%s', cell_savename);
            end
            
            if strcmp(metrics.name,'BPS_normsubtract_smooth10msec')
                eval(sprintf('load %s/%s.mat', d_load_WN,cell_savename));
                xval_WN = fittedGLM.xvalperformance; clear fittedGLM
                eval(sprintf('load %s/%s.mat', d_load_NSEM,cell_savename));
                xval_NSEM = fittedGLM.xvalperformance; clear fittedGLM
            elseif strcmp(metrics.name,'Viktor_normsubtract_25msec')
                eval(sprintf('load %s/%s.mat', d_load_WN,cell_savename));
                vksp_WN = crossval_perf; clear crossval_perf
                eval(sprintf('load %s/%s.mat', d_load_NSEM,cell_savename));
                vksp_NSEM = crossval_perf; clear crossval_perf
                
           
            end
                
             
            
            if strcmp(metrics.name,'BPS_normsubtract_smooth10msec')
                score_WN   =    xval_WN.logprob_glm_bpspike - raster_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values(rast_normindex);
                score_NSEM =  xval_NSEM.logprob_glm_bpspike - raster_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values(rast_normindex);
                
                if i_celltype == 1
                    model_comparison.scores{i_exp}.ONP_scores(1,i_cell)  = score_WN;
                    model_comparison.scores{i_exp}.ONP_scores(2,i_cell)  = score_NSEM;
                elseif i_celltype == 2
                    model_comparison.scores{i_exp}.OFFP_scores(1,i_cell) = score_WN;
                    model_comparison.scores{i_exp}.OFFP_scores(2,i_cell) = score_NSEM;
                end
            end  
            
            if strcmp(metrics.name,'Viktor_normsubtract_25msec')
                WN_timebin = find(vksp_WN.scores.Viktor_Time_Bins == 32);
                NSEM_timebin = find(vksp_NSEM.scores.Viktor_Time_Bins == 32);
                
                score_WN   =  vksp_WN.scores.metric_raw(WN_timebin) - raster_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values(rast_normindex);
                score_NSEM =  vksp_NSEM.scores.metric_raw(NSEM_timebin) - raster_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values(rast_normindex);
                
                if i_celltype == 1
                    model_comparison.scores{i_exp}.ONP_scores(1,i_cell)  = score_WN;
                    model_comparison.scores{i_exp}.ONP_scores(2,i_cell)  = score_NSEM;
                elseif i_celltype == 2
                    model_comparison.scores{i_exp}.OFFP_scores(1,i_cell) = score_WN;
                    model_comparison.scores{i_exp}.OFFP_scores(2,i_cell) = score_NSEM;
                end
            end 
        end
        
    end
    
end




%% Plotting WNWN and NSEMNSEM New

% Square Plots
minval = metrics.plots.init_minval;
maxval = metrics.plots.init_maxval;


MS = 40;
MS_ON = 20;
for i_plot = [1]
	
    for i_label = 1
        if i_label == 1, label_points = false; end
        if i_label == 2, label_points = true; end
    
    minval = Inf;
    maxval = -Inf;
    clf; hold on; set(gca,'fontsize',12);% subplot(3,3,[5,6,8,9]); 
    if i_plot == 0
        transform = 'ID'; lineplot = true;
	elseif i_plot == 1
        transform = 'ID_truncate'; lineplot = true;
    elseif i_plot == 2
        transform = 'logistic_2'; log_k=2; 
    elseif i_plot == 3
        transform = 'logistic_3'; log_k=3;
    
    end
    
    for i_exp = exptests
        if i_exp == 1, colorstring = 'r'; end
        if i_exp == 2, colorstring = 'g'; end
        if i_exp == 3, colorstring = 'b'; end
        if i_exp == 4, colorstring = 'c'; end

        for i_celltype = [2,1]

            if i_celltype == 1, vals = model_comparison.scores{i_exp}.ONP_scores; cids = model_comparison.scores{i_exp}.ONP_cells; end
            if i_celltype == 2, vals = model_comparison.scores{i_exp}.OFFP_scores; cids = model_comparison.scores{i_exp}.OFFP_cells; end

            plotstring = sprintf('%s.',colorstring);
            
            if strcmp(transform,'ID')
                vals = vals;
                vals = vals;
            end
            if length(transform) >=8 &&  strcmp(transform(1:8),'logistic')
                vals = 1./ ( 1 + exp(-(log_k)*vals) ) - .5;
                vals = 1./ ( 1 + exp(-(log_k)*vals) ) - .5;
            end
            
            if i_metric == 1 && strcmp(transform, 'ID_truncate')
                vals(find(vals<=-.75)) = -.75;
                vals(find(vals>= .25)) = .25;
            end
            
            maxval = max( maxval, max(vals(:)) );
            minval = min( minval, min(vals(:)) );           
            plot(vals(1,:), vals(2,:), plotstring, 'markersize', MS);
            if i_celltype == 1, plot(vals(1,:), vals(2,:), 'w.', 'markersize', MS_ON); end
         %   if i_celltype == 2, plot(vals(1,:), vals(2,:), 'k.', 'markersize', MS_ON); end
            if exist('label_points','var') && label_points
                for i_cell = 1:length(cids)
                    lab = sprintf('%d', cids(i_cell));
                    text(vals(1,i_cell),vals(2,i_cell),lab);
                end
            end
        end
    end
    unity_line = linspace(minval, maxval,100);
    plot(unity_line,unity_line,'k')
    xlim([minval,maxval]);
    ylim([minval,maxval]);
    if i_metric == 1 && strcmp(transform, 'ID_truncate')
        xlim([-.75,.25]);
        ylim([-.75,.25]);
        unity_line = linspace(-.75, .25,100);
        plot(unity_line,unity_line,'k')
    end
    xlabel(metrics.plots.xlabel);
    ylabel(metrics.plots.ylabel);
    
    if i_metric == 1 && strcmp(transform, 'ID_truncate')
        set(gca,'xtick', [-.75,-.5,-.25,0,.25])
        set(gca,'ytick', [-.75,-.5,-.25,0,.25])
    end
    
    plotname = sprintf('%s', transform);
    title(sprintf('%s', shortname), 'interpreter','none')

    if exist('label_points','var') && label_points
        plotname = sprintf('%s_Labels',plotname);
    end
   % orient tall
    eval(sprintf('print -dpdf %s/%s_%s.pdf', savedir, savename, plotname));
    eval(sprintf('print -deps %s/%s_%s.eps', savedir, savename, plotname));
    end
end


end
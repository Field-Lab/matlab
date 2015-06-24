% GOAL:
% Modernization of glm_AH/newtestcode/comparemodels.m
%
% CALLS:
% GLM_settings, GLM_fitname NSEM_BaseDirectories cell_list
% Standard Bookkeeping

% Version_2  Works! Automatically plots 2 versions (detailed, summary)
%            over all experiments and single populations
% Version_1  Cleaner Calling.. automatically plot everything
% Version_0  Works, produces a single plot.  2015-06-18
%  Continuation of hack_comparemodels_2


% Calling Sequences
%{
clear ; close all; clc;
%comparison_name = 'WNvsNSEM-standardGLM-noPS-LogisticfixMU'; 
%metrics = [1 2 3 4 5 6]; 
comparison_name = 'WNvsNSEM-standardGLM'; 
metrics = [1 3 4 6];
cellselection_type = 'glmconv4pct';
for i_metric = metrics
    if i_metric == 1, metric = 'BPS_divideCRM'; end
    if i_metric == 2, metric = 'VSPKD50msec_normdivide'; end
    if i_metric == 3, metric = 'FracVar10msec_normdivide'; end
    if i_metric == 4, metric = 'FracVar10msec'; end
    if i_metric == 5, metric = 'VSPKD50msec_normsubtract'; end    
    if i_metric == 6, metric = 'BPS_divideUOP'; end
    compareModels(comparison_name, metric,cellselection_type)
end


%cellselection_type = 'glmconv4pct';
%cellselection_type = 'all';
%cellselection_type = 'shortlist';

comparison_name = 'NSEM-standardGLM-PSConstrain-Sub1';



%}


function [model_comparison, outputnotes] = compareModels(comparison_name,metric,cellselection_type)
%% UNPACKING Comparison and Metic
outputnotes.read = 'Readout of which model fits still need to have scores aggregated';
outputnotes.problem_counter = 0;
exps = [1 2 3 4];
celltypes = [1 2];
plotparams.exps = exps;
plotparams.celltypes = celltypes;

% Unpack comparison_name
if strcmp(comparison_name, 'WNvsNSEM-standardGLM')
    models{1}.settings = {};
    models{1}.fit_type = 'WN';
    models{2}.settings = {};
    models{2}.fit_type = 'NSEM';
    plotparams.xlabel             = 'White Noise';
    plotparams.ylabel             = 'Natural Scenes';
    plotparams.title_comparison   = 'Standard GLM (no coupling)';
    plotparams.purpose            = 'Check WN vs NSEM for standard GLM (fixed space to WN-STA, fit time, fit post-spike, no coupling)';
end
if strcmp(comparison_name, 'WNvsNSEM-standardGLM-noPS-LogisticfixMU')
    models{1}.settings{1}.type = 'PostSpikeFilter';
    models{1}.settings{1}.name =  'OFF';
    models{1}.fit_type = 'WN';
    models{1}.postfilerNL = 'Logistic_fixMU';
    models{2}.settings{1}.type = 'PostSpikeFilter';
    models{2}.settings{1}.name =  'OFF';
    models{2}.fit_type = 'NSEM';
    models{2}.postfilerNL         = 'Logistic_fixMU';
    plotparams.xlabel             = 'White Noise';
    plotparams.ylabel             = 'Natural Scenes';
    plotparams.title_comparison   = 'LN with Logistic';
    plotparams.purpose            = 'Check WN vs NSEM for a general LN model (GLM no PS then fit logistic)';
end

if strcmp(comparison_name, 'NSEM-standardGLM-PSConstrain-Sub1')
    models{1}.settings = {};
    models{1}.fit_type = 'NSEM';
    models{2}.settings= {};
    models{2}.fit_type = 'NSEM';
    models{2}.postfilerNL         = 'PS_Contrain_sub1';
    plotparams.xlabel             = 'White Noise';
    plotparams.ylabel             = 'Natural Scenes';
    plotparams.title_comparison   = 'LN with Logistic';
    plotparams.purpose            = 'Check WN vs NSEM for a general LN model (GLM no PS then fit logistic)';
end


% Unpack metric
switch metric
    case 'BPS_divideUOP'
        rawmetric_name = 'crossval_BPS';
        normalize.doit = true;
        normalize.name = 'base_crm_findPS';
        normalize.extension = 'uop_bps';
        normalize.operation = 'divide';        
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        plotparams.title_metric = 'Fractional Bits';
        
    case 'BPS_divideCRM'
        rawmetric_name = 'crossval_BPS';
        normalize.doit = true;
        normalize.name = 'base_crm_findPS';
        normalize.extension = 'crm_bps';
        normalize.operation = 'divide';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        plotparams.title_metric = 'Fractional Bits';
        
    case 'FracVar10msec'
        rawmetric_name = 'crossval_fracvar_10msec';
        normalize.doit = false;
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        plotparams.title_metric = 'Fraction of Variance';
        
    case 'FracVar10msec_normdivide'
        rawmetric_name = 'crossval_fracvar_10msec';
        normalize.doit = true;
        normalize.name = 'fracvar_10msec_ODDEVEN';
        normalize.extension = '';
        normalize.operation = 'divide';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        plotparams.title_metric = 'Normed Fraction of Variance';
        
    case 'VSPKD50msec'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = false;
        plotparams.low_lim = 0;
        plotparams.high_lim = 2;
        plotparams.title_metric = 'Victor Spike (50msec)';
  
    case 'VSPKD50msec_normsubtract'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = true;
        normalize.name = 'vspkd_msec50';
        normalize.extension = 'rast_vspkd';
        normalize.operation = 'subtract';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        plotparams.title_metric = 'NormSubtract Victor Spike (50msec)';
        
    case 'VSPKD50msec_normdivide'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = true;
        normalize.name = 'vspkd_msec50';
        normalize.extension = 'rast_vspkd';
        normalize.operation = 'divide';
        plotparams.low_lim = 1;
        plotparams.high_lim = 2;
        plotparams.title_metric = 'Normalized Victor Spike (50msec)';
end
plotparams.title = sprintf('%s: %s', plotparams.title_metric, plotparams.title_comparison);
%% Bookkeeping
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
if strcmp(cellselection_type, 'glmconv4pct')
    eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
end
currentdir = pwd;
savedir = sprintf('%s/Plots/ModelComparisons/%s/%s',BD.GLM_output_analysis,comparison_name,metric);
if ~exist(savedir,'dir'), mkdir(savedir); end

savename = sprintf('%s_%s_%s',comparison_name,metric,cellselection_type);
clear expstring celltypestring

models{1}.GLMType         = GLM_settings('default',models{1}.settings);
models{1}.fitname         = GLM_fitname(models{1}.GLMType);
if isfield(models{1}, 'postfilterNL') && strcmp(models{1}.postfilerNL,'Logistic_fixMU')
    models{1}.fitname = sprintf('%s/Logistic_fixMU', models{1}.fitname);
end
models{2}.GLMType         = GLM_settings('default',models{2}.settings);
models{2}.fitname         = GLM_fitname(models{2}.GLMType);
if isfield(models{1}, 'postfilterNL') && strcmp(models{2}.postfilerNL,'Logistic_fixMU')
    models{2}.fitname = sprintf('%s/Logistic_fixMU', models{2}.fitname);
end

models{1}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, models{1}.fitname);
models{2}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, models{2}.fitname);

%%% LOADING 
model_comparison.byexpnm         = allcells;
model_comparison.comparison_name = comparison_name;
model_comparison.metric          = metric;
model_comparison.models          = models;
model_comparison.timestamp       = datestr(clock);
model_comparison.code_name       = mfilename('fullpath');
%% GET CELLS SQUARED UP
% Define Cell_Subset structure
for i_exp = exps
    cell_subset = cell(2,1);
    cell_subset{1}.celltype = 'ONPar';
    cell_subset{1}.all_cids = allcells{i_exp}.ONP;
    cell_subset{2}.celltype = 'OFFPar';
    cell_subset{2}.all_cids = allcells{i_exp}.OFFP;

    % Grab cids of cell subsets
    if strcmp(cellselection_type, 'shortlist')
        [~,candidate_cells,~]  = cell_list(i_exp,cellselection_type);
        candidate_cells = cell2mat(candidate_cells);
        for i_celltype = 1:2
            cell_subset{i_celltype}.subset_cids = intersect(candidate_cells, cell_subset{i_celltype}.all_cids);
        end
        clear candidate_cells
    elseif strcmp(cellselection_type,'glmconv4pct')
        conv_column     = 2; 
        conv_index_ON   = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
        conv_index_OFF  = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
        cell_subset{1}.subset_cids = allcells{i_exp}.ONP(conv_index_ON);
        cell_subset{2}.subset_cids = allcells{i_exp}.OFFP(conv_index_OFF);
        clear conv_index_ON conv_index_OFF conv_columns
    elseif strcmp(cellselection_type,'all')
        cell_subset{1}.subset_cids = allcells{i_exp}.ONP;
        cell_subset{2}.subset_cids = allcells{i_exp}.OFFP;
    end
    
    % Translate into indices
    for i_celltype = 1:2
        full   = cell_subset{i_celltype}.all_cids;
        subset = cell_subset{i_celltype}.subset_cids;
        index  = NaN(size(subset));
        for i_cell = 1:length(subset)
            cid = subset(i_cell);
            index(i_cell) = find(full == cid);
        end
        cell_subset{i_celltype}.subset_indices = index;
        clear index subset cid full
    end
    model_comparison.byexpnm{i_exp}.cell_subset =cell_subset; 
end
%% Build up data
for i_exp = exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    
    % Pull out relevant structures for tracking results
    cell_subset       = model_comparison.byexpnm{i_exp}.cell_subset;
    scores_bycelltype = cell_subset;
    
    % Update with each new metric here!
    for i_model = 1:2
        clear aggregated_scores raster_scores
        
        scores_filename = sprintf('%s/%s_%s.mat',models{i_model}.aggscores_dir, rawmetric_name,exp_nm); 
        if ~exist(scores_filename)
            outputnotes.problem_counter = outputnotes.problem_counter+1;
            message = sprintf('Aggregated Scores Structure is missing.  Experiment: %s, Model: %s',...
                exp_nm, models{i_model}.fitname);
            outputnotes.notes{outputnotes.problem_counter} = message;
            display(message); clear message
            plotparams.exps = setdiff(plotparams.exps,i_exp);
            continue
        end
            
        eval(sprintf('load %s/%s_%s.mat aggregated_scores',  models{i_model}.aggscores_dir, rawmetric_name,exp_nm));
        if normalize.doit
            eval(sprintf('load %s/Raster_Metrics/%s_%s.mat raster_scores', BD.Cell_Selection, normalize.name, exp_nm));
        end

        for i_celltype = celltypes
            scores_bycelltype{i_celltype}.cids = cell_subset{i_celltype}.subset_cids;
            
            % no way around it
            if strcmp(models{i_model}.fit_type, 'WN')
                rawscores_all = aggregated_scores.celltype{i_celltype}.scores_WN;
            elseif strcmp(models{i_model}.fit_type, 'NSEM')
                rawscores_all = aggregated_scores.celltype{i_celltype}.scores_NSEM;
            end
            rawscores_subset = rawscores_all(cell_subset{i_celltype}.subset_indices);
            
            if normalize.doit
                if strcmp(models{i_model}.fit_type, 'WN')
                    normalizers_struct = raster_scores.celltype{i_celltype}.scores_WN;
                elseif strcmp(models{i_model}.fit_type, 'NSEM')
                    normalizers_struct = raster_scores.celltype{i_celltype}.scores_NSEM;
                end
                % super hack no way around it for now 
                if strcmp(normalize.extension,'uop_bps')
                    normalizers_all = normalizers_struct.uop_bps;
                elseif strcmp(normalize.extension,'crm_bps')
                    normalizers_all = normalizers_struct.crm_bps;
                elseif strcmp(normalize.extension,'crm_bps')
                    normalizers_all = normalizers_struct.crm_bps;
                elseif strcmp(normalize.extension,'rast_vspkd')
                    normalizers_all = normalizers_struct.rast_vspkd;
                elseif isempty(normalize.extension)
                    normalizers_all = normalizers_struct;
                else
                    error('Normalizers Extension is not Properly Defined!');
                end
                normalizers_subset = normalizers_all(cell_subset{i_celltype}.subset_indices);
            end
           
            % Define the final score
            if normalize.doit
                switch normalize.operation
                    case 'divide'
                        finalscores = rawscores_subset./(normalizers_subset); 
                    case 'subtract'
                        finalscores = rawscores_subset - (normalizers_subset);                    
                end
            else
                finalscores = rawscores_subset;
            end
            
            % Safety:  Track Wierd Numbers
            if sum(isnan(finalscores) + isinf(finalscores)) > 0
                problem_indices = find(isnan(finalscores) + isinf(finalscores));
                cell_list_str = 'Cells';
                for i_ind = 1:length(problem_indices)
                    p_ind = problem_indices(i_ind);
                    p_cid = cell_subset{i_celltype}.subset_cids(p_ind);
                    cell_list_str = sprintf('%s, %s', cell_list_str, num2str(p_cid));
                end
                clear p_ind p_cid
                
                if i_celltype == 1 
                    frontmessage = sprintf('Metric is NAN or INF for Exp: %s, ONParsol:', exp_nm);
                elseif i_celltype == 2
                    frontmessage = sprintf('Metric is NAN or INF for Exp: %s, OFFParsols:',exp_nm);
                end
                
                outputnotes.problem_counter = outputnotes.problem_counter+1;
                message = sprintf('%s %s ::: Model: %s ',...
                    frontmessage, cell_list_str, models{i_model}.fitname);
                outputnotes.notes{outputnotes.problem_counter} = message;
                display(message); clear message
                clear frontmessage cell_list_str
                
                plotnotes = sprintf('Problem Scores exist in experiment %d', i_exp)
            end
            
            
            if i_model == 1, scores_bycelltype{i_celltype}.finalscores_model1 = finalscores; end
            if i_model == 2, scores_bycelltype{i_celltype}.finalscores_model2 = finalscores; end           
        end
    end
    
    % Put scores structure back in
    model_comparison.byexpnm{i_exp}.scores_bycelltype =  scores_bycelltype;
end
eval(sprintf('save %s/data_%s.mat model_comparison outputnotes', savedir, savename))

%% plotting section

expstring = 'exps';
if find(plotparams.exps ==1), expstring = sprintf('%sA',expstring); end
if find(plotparams.exps ==2), expstring = sprintf('%sB',expstring); end
if find(plotparams.exps ==3), expstring = sprintf('%sC',expstring); end
if find(plotparams.exps ==4), expstring = sprintf('%sD',expstring); end
if strcmp(expstring, 'expsABCD')
    plotparams.exps = [2 3 1 4];
end
printname_base = sprintf('%s_%s',savename,expstring);
% hack to prevent pdf printing error
homedir = pwd;

printname_notext = sprintf('%s_NOTEXT',printname_base);
printname_fullplot = sprintf('%s_FULLPLOTS',printname_base);

cd(savedir)
subR_plotcomparison_nolabels(model_comparison,plotparams,printname_notext);
subR_plotcomparison_fullplots(model_comparison,plotparams,printname_fullplot);
cd(homedir)

% hack to get plot without the 2013-10-10-0
if strcmp(expstring, 'expsABCD')
    expstring_no4  = 'expsABC';
    plotparams_no4 = plotparams;
    plotparams_no4.exps = [2 3 1];
    printname_notext   = sprintf('%s_%s_NOTEXT',savename,expstring_no4);
    printname_fullplot = sprintf('%s_%s_FULLPLOTS',savename,expstring_no4);
    
    cd(savedir)
    subR_plotcomparison_nolabels(model_comparison,plotparams_no4,printname_notext);
    subR_plotcomparison_fullplots(model_comparison,plotparams_no4,printname_fullplot);
    cd(homedir)
    clear expstring_no4 plotparams_no4
end

subset_params = plotparams;
for i_exp = 1:length(plotparams.exps)
    subset_params.exps = plotparams.exps(i_exp);
    expstring = 'exps';
    if find(subset_params.exps ==1), expstring = sprintf('%sA',expstring); end
    if find(subset_params.exps ==2), expstring = sprintf('%sB',expstring); end
    if find(subset_params.exps ==3), expstring = sprintf('%sC',expstring); end
    if find(subset_params.exps ==4), expstring = sprintf('%sD',expstring); end

    for i_celltype = 1:2
        subset_params.celltypes = i_celltype;
        if i_celltype == 1, celltypestring = 'ONParasols'; end
        if i_celltype == 2, celltypestring = 'OFFParasols'; end
        printname_notext = sprintf('%s_%s_%s_NOTEXT',savename,expstring,celltypestring);
        printname_fullplot = sprintf('%s_%s_%s_FULLPLOTS',savename,expstring,celltypestring);
        cd(savedir)
        subR_plotcomparison_nolabels(model_comparison,subset_params,printname_notext);
        subR_plotcomparison_fullplots(model_comparison,subset_params,printname_fullplot);
        cd(homedir)
        
    end
end

end


function subR_plotcomparison_nolabels(model_comparison,plotparams,printname)
% subRoutine form done 2015-06-18 AKHeitman
% PreSet plotting limits, clean colors etc.
% Main Call: only reference to model_comparison structure 
% scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
% scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;
% plot_entries

colors = {'r.','g.','b.','c.'};
MS_A = 12;
MS_B = 30; 
clf; hold on
axis square
low_lim  = -Inf;
high_lim = Inf;

if isfield(plotparams, 'low_lim')
    low_lim = plotparams.low_lim;
end
if isfield(plotparams, 'high_lim')
    high_lim = plotparams.high_lim;
end

xlim([low_lim, high_lim]);
ylim([low_lim, high_lim]);
plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');

for i_exp = plotparams.exps
    colorstring = colors{i_exp};
    for i_celltype = plotparams.celltypes
        scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
        scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;

        scores1(find(scores1<=low_lim)) = low_lim;
        scores2(find(scores2<=low_lim)) = low_lim;

        scores1(find(scores1>=high_lim)) = high_lim;
        scores2(find(scores2>=high_lim)) = high_lim;
        
        plot(scores1, scores2, colorstring, 'markersize', MS_B);
        if i_celltype == 1
            plot(scores1, scores2, 'w.', 'markersize', MS_A);
        elseif i_celltype == 2
            plot(scores1, scores2, 'k.', 'markersize', MS_A);
        end            
    end
end
orient portrait
eval(sprintf('print -dpdf plot_%s.pdf',printname))
end


function subR_plotcomparison_fullplots(model_comparison,plotparams,printname)
% subRoutine form done 2015-06-18 AKHeitman
% PreSet plotting limits, clean colors etc.
% Main Call: only reference to model_comparison structure 
% scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
% scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;
% plot_entries

clf
subplot(3,1,1); axis off
delta = .15;
c=0;   text(-.1, 1-delta*c,sprintf('PURPOSE: %s', plotparams.purpose ));
c=c+1; text(-.1, 1-delta*c,sprintf('Metric: %s,  Comparison Title: %s',...
    model_comparison.metric, model_comparison.comparison_name ),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('X-axis: %s Fit: %s',...
    model_comparison.models{1}.fit_type,model_comparison.models{1}.fitname),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('Y-axis: %s Fit: %s',...
    model_comparison.models{2}.fit_type,model_comparison.models{2}.fitname),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,...
    'Outer Dot Color: { (Red: 2012-08-09-3), (Green: 2012-09-27-3), (Blue: 2013-08-19-6), (Cyan: 2013-10-10-0) }','interpreter','none');
c=c+1; text(-.1, 1-delta*c,...
    'Inner Dot Color: { (White: On Parasol), (Black: Off Parasol) }','interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-delta*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none' );
    


colors = {'r.','g.','b.','c.'};
MS_A = 10;
MS_B = 24;
low_lim  = -Inf;
high_lim = Inf;
if isfield(plotparams, 'low_lim')
    low_lim = plotparams.low_lim;
end
if isfield(plotparams, 'high_lim')
    high_lim = plotparams.high_lim;
end





subplot(3,2,[3 5]);
hold on
axis square
minval = Inf; maxval = -Inf;
set(gca, 'fontsize', 10);
xlabel(plotparams.xlabel);
ylabel(plotparams.ylabel);
title(sprintf('%s: Raw Plot', plotparams.title)); 
for i_exp = plotparams.exps
    for i_celltype = plotparams.celltypes
        scores = [model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1; ...
        model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2];
        minval = min(minval, min(scores(:)) );
        maxval = max(maxval, max(scores(:)) );
    end
end
plot(linspace(minval,maxval,100),linspace(minval,maxval,100),'k');
xlim([minval, maxval]);
ylim([minval, maxval]);
for i_exp = plotparams.exps
    colorstring = colors{i_exp};
    for i_celltype = plotparams.celltypes
        scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
        scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;
        plot(scores1, scores2, colorstring, 'markersize', MS_B);
        if i_celltype == 1
            plot(scores1, scores2, 'w.', 'markersize', MS_A);
        elseif i_celltype == 2
            plot(scores1, scores2, 'k.', 'markersize', MS_A);
        end            
    end
end



subplot(3,2,[4 6]);
hold on
axis square
set(gca, 'fontsize', 10);
xlim([low_lim, high_lim]);
ylim([low_lim, high_lim]);
xlabel(plotparams.xlabel);
ylabel(plotparams.ylabel);
title(sprintf('%s: Truncate Extreme Values', plotparams.title)); 
plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');
for i_exp = plotparams.exps
    colorstring = colors{i_exp};
    for i_celltype = plotparams.celltypes
        scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
        scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;

        scores1(find(scores1<=low_lim)) = low_lim;
        scores2(find(scores2<=low_lim)) = low_lim;

        scores1(find(scores1>=high_lim)) = high_lim;
        scores2(find(scores2>=high_lim)) = high_lim;
        
        plot(scores1, scores2, colorstring, 'markersize', MS_B);
        if i_celltype == 1
            plot(scores1, scores2, 'w.', 'markersize', MS_A);
        elseif i_celltype == 2
            plot(scores1, scores2, 'k.', 'markersize', MS_A);
        end            
    end
end

orient landscape
eval(sprintf('print -dpdf plot_%s.pdf',printname))
end
  



  


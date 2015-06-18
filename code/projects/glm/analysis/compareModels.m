% Version 2 last hacked version that works cleanly
% Version 1 use the NEW STATS
% Version 0 only looked at fittedGLM and not other metrics 2015-03-18
% Version 0 more compatible with shortlist types
% hack_comparemodels
% AKHEITMAN 2015-03-13
% 
% GOAL:
% Modernization of glm_AH/newtestcode/comparemodels.m
%
% CALLS:
% GLM_settings, GLM_fitname NSEM_BaseDirectories cell_list
% Standard Bookkeeping
clear ; close all; clc;

% Dictate comparison plot
%{
comparison_name = 'WNvsNSEM_standardGLM_noPS';
details.metric = 'BPS_divideUOP';
%details.metric = 'BPS_divideCRM';
%details.metric = 'FracVar_10msec'
%details.metric = 'FracVar_20msec';
%details.metric = 'VSPKD_50msec';
%details.metric = 'VSPKD_50msec_subtract';
%details.metric = 'VSPKD_50msec_divide';
model{1}.settings{1}.type = 'PostSpikeFilter';
model{1}.settings{1}.name =  'OFF';
model{1}.fit_type = 'WN';

model{1}.special_arg = '';
model{2}.settings{1}.type = 'PostSpikeFilter';
model{2}.settings{1}.name =  'OFF';
model{2}.fit_type = 'NSEM';
model{2}.special_arg = '';
%}


comparison_name = 'WNvsNSEM_standardGLM_noPS_LOG';
%details.metric = 'BPS_divideUOP';
%details.metric = 'BPS_divideCRM';
details.metric = 'FracVar_10msec'
%details.metric = 'FracVar_20msec';
%details.metric = 'VSPKD_50msec';
%details.metric = 'VSPKD_50msec_subtract';
%details.metric = 'VSPKD_50msec_divide';

%{
comparison_name = 'WNvsNSEM_standardGLM';
model{1}.settings= {};
model{1}.fit_type = 'WN';
model{2}.settings= {};
model{2}.fit_type = 'NSEM';
details.metric = 'BPS_divideCRM'
%}


%details.cellselection_type = 'all';
details.cellselection_type = 'glmconv_4pct';
details.celltypes = [1 2];
details.exps = [4];

function compareModels(comparison_name, metric, cellselection_type)

% Unpack comparison_name
if strcmp(comparison_name, 'WNvsNSEM_standardGLM_noPS_LOG')
    models{1}.settings{1}.type = 'PostSpikeFilter';
    models{1}.settings{1}.name =  'OFF';
    models{1}.fit_type = 'WN';
    models{1}.special_arg = 'Logistic_fixMU';
    models{2}.settings{1}.type = 'PostSpikeFilter';
    models{2}.settings{1}.name =  'OFF';
    models{2}.fit_type = 'NSEM';
    models{2}.special_arg = 'Logistic_fixMU';
end



%% Bookkeeping
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
if strcmp(details.cellselection_type, 'glmconv_4pct')
    eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
end
currentdir = pwd;
savedir = sprintf('%s/Plots/ModelComparisons/%s',BD.GLM_output_analysis,comparison_name);
if ~exist(savedir,'dir'), mkdir(savedir); end

expstring = 'exp';
if find(details.exps ==1), expstring = sprintf('%sA',expstring); end
if find(details.exps ==2), expstring = sprintf('%sB',expstring); end
if find(details.exps ==3), expstring = sprintf('%sC',expstring); end
if find(details.exps ==4), expstring = sprintf('%sD',expstring); end
celltypestring = '';
if length(details.celltypes) == 1
    if details.celltypes == 1, celltypestring = 'ONonly'; end
    if details.celltypes == 2, celltypestring = 'OFFonly'; end
end
savename = sprintf('%s_%s_%s_%s%s',comparison_name,details.metric, details.cellselection_type, expstring, celltypestring);
clear expstring celltypestring

model{1}.GLMType         = GLM_settings('default',model{1}.settings);
model{1}.fitname         = GLM_fitname(model{1}.GLMType);

if strcmp(model{1}.special_arg,'Logistic_fixMU')
    model{1}.fitname = sprintf('%s/Logistic_fixMU', model{1}.fitname);
end


model{1}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, model{1}.fitname);
model{2}.GLMType         = GLM_settings('default',model{2}.settings);
model{2}.fitname         = GLM_fitname(model{2}.GLMType);
model{2}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, model{2}.fitname);

if strcmp(model{2}.special_arg,'Logistic_fixMU')
    model{2}.fitname = sprintf('%s/Logistic_fixMU', model{2}.fitname);
end

%%% LOADING 
model_comparison.byexpnm         = allcells;
model_comparison.comparison_name = comparison_name;
model_comparison.notes.n1        = 'Use Normalization from ';
model_comparison.notes.n2        = 'UOP: unconditioned optimal performance';
model_comparison.notes.n3        = 'CRM: conditioned rate model';
model_comparison.notes.n4        = 'BPS: Bits Per Spike';
model_comparison.notes.n5        = 'LPPS: logarithmic probability per second';
model_comparison.notes.n6        = 'GLM: Generalized Linear Model';
model_comparison.timestamp       = datestr(clock);
model_comparison.code_name       = mfilename('fullpath');
%% GET CELLS SQUARED UP
% Define Cell_Subset structure
for i_exp = details.exps
    cell_subset = cell(2,1);
    cell_subset{1}.celltype = 'ONPar';
    cell_subset{1}.all_cids = allcells{i_exp}.ONP;
    cell_subset{2}.celltype = 'OFFPar';
    cell_subset{2}.all_cids = allcells{i_exp}.OFFP;

    % Grab cids of cell subsets
    if strcmp(details.cellselection_type, 'shortlist')
        [~,candidate_cells,~]  = cell_list(i_exp,cellselection_type);
        candidate_cells = cell2mat(candidate_cells);
        for i_celltype = 1:2
            cell_subset{i_celltype}.subset_cids = intersect(candidate_cells, cell_subset{i_celltype}.all_cids);
        end
        clear candidate_cells
    elseif strcmp(details.cellselection_type,'glmconv_4pct')
        conv_column     = 2; 
        conv_index_ON   = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
        conv_index_OFF  = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
        cell_subset{1}.subset_cids = allcells{i_exp}.ONP(conv_index_ON);
        cell_subset{2}.subset_cids = allcells{i_exp}.OFFP(conv_index_OFF);
        clear conv_index_ON conv_index_OFF conv_columns
    elseif strcmp(details.cellselection_type,'all')
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
for i_exp = details.exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    % Load normalizer by experiment number
    if strcmp(details.metric, 'BPS_divideUOP') || strcmp(details.metric, 'BPS_divideCRM')
        neednormalization = true;
        if i_exp < 4
            normalizer = 'base_crm_findPS';
        elseif i_exp == 4
            normalizer = 'base_crm_importPS';
        end    
        eval(sprintf('load %s/Raster_Metrics/%s_%s.mat raster_scores', BD.Cell_Selection, normalizer, exp_nm));
        clear normalizer
    elseif strcmp(details.metric,'VSPKD_50msec_subtract') || strcmp(details.metric,'VSPKD_50msec_divide') 
        neednormalization = true;
        normalizer = 'vspkd_msec50'
        eval(sprintf('load  %s/Raster_Metrics/%s_%s.mat raster_scores', BD.Cell_Selection, normalizer, exp_nm));
        clear normalizer
    else
        neednormalization = false;
    end
    
    % Pull out relevant structures for tracking results
    cell_subset       = model_comparison.byexpnm{i_exp}.cell_subset;
    scores_bycelltype = cell_subset;
    
    
    % Update with each new metric here!
    for i_model = 1:2
        if strcmp(details.metric,'BPS_divideUOP'),  underlyingmetric_name = sprintf('crossval_BPS_%s',exp_nm); end
        if strcmp(details.metric,'BPS_divideCRM'),  underlyingmetric_name = sprintf('crossval_BPS_%s',exp_nm); end  
        if strcmp(details.metric,'FracVar_10msec'), underlyingmetric_name = sprintf('crossval_fracvar_10msec_%s',exp_nm); end
        if strcmp(details.metric,'FracVar_20msec'), underlyingmetric_name = sprintf('crossval_fracvar_20msec_%s',exp_nm); end
        if strcmp(details.metric,'VSPKD_50msec'),   underlyingmetric_name = sprintf('crossval_victorspike_50msec_%s',exp_nm); end
        if strcmp(details.metric,'VSPKD_50msec_subtract'), underlyingmetric_name = sprintf('crossval_victorspike_50msec_%s',exp_nm); end
        if strcmp(details.metric,'VSPKD_50msec_divide'), underlyingmetric_name = sprintf('crossval_victorspike_50msec_%s',exp_nm); end
        
        eval(sprintf('load %s/%s.mat aggregated_scores',  model{i_model}.aggscores_dir, underlyingmetric_name));
        
        for i_celltype = details.celltypes
            scores_bycelltype{i_celltype}.cids = cell_subset{i_celltype}.subset_cids;
            if strcmp(model{i_model}.fit_type, 'WN')
                rawscores_all = aggregated_scores.celltype{i_celltype}.scores_WN;
                if neednormalization
                    normalizers_all = raster_scores.celltype{i_celltype}.scores_WN;
                end
            elseif strcmp(model{i_model}.fit_type, 'NSEM')
                rawscores_all = aggregated_scores.celltype{i_celltype}.scores_NSEM;
                if neednormalization
                    normalizers_all = raster_scores.celltype{i_celltype}.scores_NSEM;
                end
            end
            rawscores_subset = rawscores_all(cell_subset{i_celltype}.subset_indices);
            
            % Define the final score
            if strcmp(details.metric,'BPS_divideUOP')
                finalscores = rawscores_subset./(normalizers_all.uop_bps(cell_subset{i_celltype}.subset_indices)); 
            elseif strcmp(details.metric, 'BPS_divideCRM')
                finalscores = rawscores_subset./(normalizers_all.crm_bps(cell_subset{i_celltype}.subset_indices));
            elseif strcmp(details.metric,'FracVar_10msec')
                finalscores = rawscores_subset;
            elseif strcmp(details.metric,'FracVar_20msec')
                finalscores = rawscores_subset;
            elseif strcmp(details.metric,'VSPKD_50msec')
                finalscores = rawscores_subset;
            elseif strcmp(details.metric,'VSPKD_50msec_subtract')
                finalscores = rawscores_subset - (normalizers_all.rast_vspkd(cell_subset{i_celltype}.subset_indices));   
            elseif strcmp(details.metric,'VSPKD_50msec_divide')
                finalscores = rawscores_subset./(normalizers_all.rast_vspkd(cell_subset{i_celltype}.subset_indices));   
            end
            
            if i_model == 1, scores_bycelltype{i_celltype}.finalscores_model1 = finalscores; end
            if i_model == 2, scores_bycelltype{i_celltype}.finalscores_model2 = finalscores; end           
        end
    end
    
    % Put scores structure back in
    model_comparison.byexpnm{i_exp}.scores_bycelltype =  scores_bycelltype;
end

eval(sprintf('save %s/%s.mat model_comparison', savedir, savename))

%% PLOT CHANGE IN WN VS NSEM  (ONE SINGLE LARGE PLOT)
% hack to prevent long pdf name which can cause error
cd(savedir)


function subR_plot_comparison_nolabels(model_comparison,metric,savename,plot_entries)
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
if strcmp(metric,'BPS_divideCRM') || strcmp(metric,'BPS_divideUOP') || strcmp(metric, 'FracVar_10msec')
    low_lim  = 0;
    high_lim = 1;
end
if strcmp(metric, 'FracVar_10msec') || strcmp(metric, 'FracVar_20msec')
    low_lim  = 0;
    high_lim = 1;
end
if strcmp(metric,'VSPKD_50msec')
    low_lim = 0;
    high_lim = 2;
end
if strcmp(metric,'VSPKD_50msec_subtract')
    low_lim = 0;
    high_lim = 1;
end
if strcmp(metric,'VSPKD_50msec_divide')
    low_lim  = 1;
    high_lim = 2;
end
xlim([low_lim, high_lim]);
ylim([low_lim, high_lim]);
plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');

for i_exp = plot_entries.exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    colorstring = colors{i_exp};
    for i_celltype = plot_entries.celltypes
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
eval(sprintf('print -dpdf plot_%s.pdf',savename))
end
  


% AKHEITMAN 2015-03-13
% 
% GOAL:
% Modernization of glm_AH/newtestcode/comparemodels.m
%
% CALLS:
% GLM_settings, GLM_fitname NSEM_BaseDirectories cell_list
% Standard Bookkeeping


% Version_0  Works, produces a single plot.  2015-06-18
%  Continuation of hack_comparemodels_2



%{
clear ; close all; clc;
comparison_name = 'WNvsNSEM-standardGLM-noPS-LogisticfixMU';
%metric = 'BPS_divideUOP';
%metric = 'VSPKD50msec_normdivide';
%metric='VSPKD50msec_normsubtract';
%metric= 'FracVar10msec';
metric = 'FracVar10msec_normdivide';
cellselection_type = 'glmconv4pct';
compareModels(comparison_name, metric,cellselection_type)





%metric,'BPS_divideCRM') 
%metric,'BPS_divideUOP')
%metric, 'FracVar10msec')
%metric, 'FracVar20msec')
%metric,'VSPKD50msec')
%metric,'VSPKD50msec_normsubtract'
%metric,'VSPKD50msec_normdivide'

%}


function compareModels(comparison_name, metric, cellselection_type)
%% UNPACKING Comparison and Metic

% Unpack comparison_name
if strcmp(comparison_name, 'WNvsNSEM-standardGLM-noPS-LogisticfixMU')
    models{1}.settings{1}.type = 'PostSpikeFilter';
    models{1}.settings{1}.name =  'OFF';
    models{1}.fit_type = 'WN';
    models{1}.postfilerNL = 'Logistic_fixMU';
    models{2}.settings{1}.type = 'PostSpikeFilter';
    models{2}.settings{1}.name =  'OFF';
    models{2}.fit_type = 'NSEM';
    models{2}.postfilerNL = 'Logistic_fixMU';
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
        
    case 'BPS_divideCRM'
        rawmetric_name = 'crossval_BPS';
        normalize.doit = true;
        normalize.name = 'base_crm_findPS';
        normalize.extension = 'crm_bps';
        normalize.operation = 'divide';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        
    case 'FracVar10msec'
        rawmetric_name = 'crossval_fracvar_10msec';
        normalize.doit = false;
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        
    case 'FracVar10msec_normdivide'
        rawmetric_name = 'crossval_fracvar_10msec';
        normalize.doit = true;
        normalize.name = 'fracvar_10msec_ODDEVEN';
        normalize.extension = '';
        normalize.operation = 'divide';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        
    case 'VSPKD50msec'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = false;
        plotparams.low_lim = 0;
        plotparams.high_lim = 2;
  
    case 'VSPKD50msec_normsubtract'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = true;
        normalize.name = 'vspkd_msec50';
        normalize.extension = 'rast_vspkd';
        normalize.operation = 'subtract';
        plotparams.low_lim = 0;
        plotparams.high_lim = 1;
        
    case 'VSPKD50msec_normdivide'
        rawmetric_name = 'crossval_victorspike_50msec';
        normalize.doit = true;
        normalize.name = 'vspkd_msec50';
        normalize.extension = 'rast_vspkd';
        normalize.operation = 'divide';
        plotparams.low_lim = 1;
        plotparams.high_lim = 2;
end

%% Bookkeeping
exps = [1 2 3 4];
celltypes = [1 2];
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
if strcmp(cellselection_type, 'glmconv4pct')
    eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
end
currentdir = pwd;
savedir = sprintf('%s/Plots/ModelComparisons/%s',BD.GLM_output_analysis,comparison_name);
if ~exist(savedir,'dir'), mkdir(savedir); end

expstring = 'exp';
if find(exps ==1), expstring = sprintf('%sA',expstring); end
if find(exps ==2), expstring = sprintf('%sB',expstring); end
if find(exps ==3), expstring = sprintf('%sC',expstring); end
if find(exps ==4), expstring = sprintf('%sD',expstring); end
celltypestring = '';
if length(celltypes) == 1
    if celltypes == 1, celltypestring = 'ONonly'; end
    if celltypes == 2, celltypestring = 'OFFonly'; end
end
savename = sprintf('%s_%s_CELLS-%s%s%s',comparison_name,metric,cellselection_type, expstring, celltypestring);
clear expstring celltypestring

models{1}.GLMType         = GLM_settings('default',models{1}.settings);
models{1}.fitname         = GLM_fitname(models{1}.GLMType);
if strcmp(models{1}.postfilerNL,'Logistic_fixMU')
    models{1}.fitname = sprintf('%s/Logistic_fixMU', models{1}.fitname);
end
models{2}.GLMType         = GLM_settings('default',models{2}.settings);
models{2}.fitname         = GLM_fitname(models{2}.GLMType);
if strcmp(models{2}.postfilerNL,'Logistic_fixMU')
    models{2}.fitname = sprintf('%s/Logistic_fixMU', models{2}.fitname);
end

models{1}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, models{1}.fitname);
models{2}.aggscores_dir   = sprintf('%s/%s', BD.GLM_output_analysis, models{2}.fitname);

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
            if i_model == 1, scores_bycelltype{i_celltype}.finalscores_model1 = finalscores; end
            if i_model == 2, scores_bycelltype{i_celltype}.finalscores_model2 = finalscores; end           
        end
    end
    
    % Put scores structure back in
    model_comparison.byexpnm{i_exp}.scores_bycelltype =  scores_bycelltype;
end
eval(sprintf('save %s/%s.mat model_comparison', savedir, savename))

%%
thisdir = pwd;
plotparams.exps = [1 2 3 4];
plotparams.celltypes = [1 2];
savename = savename;
metric = 'BPS_divideUOP';
subR_plotcomparison_nolabels(model_comparison,plotparams,savedir,savename)
cd(thisdir)

end


function subR_plotcomparison_nolabels(model_comparison,plotparams,savedir,savename)
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

% hack to prevent pdf printing error
cd(savedir)
eval(sprintf('print -dpdf plot_%s.pdf',savename))
end
  


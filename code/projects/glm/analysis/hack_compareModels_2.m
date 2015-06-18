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
model{1}.settings{1}.type = 'PostSpikeFilter';
model{1}.settings{1}.name =  'OFF';
model{1}.fit_type = 'WN';

model{1}.special_arg = '';
model{2}.settings{1}.type = 'PostSpikeFilter';
model{2}.settings{1}.name =  'OFF';
model{2}.fit_type = 'NSEM';
model{2}.special_arg = '';
model{1}.special_arg = 'Logistic_fixMU';
model{2}.special_arg = 'Logistic_fixMU';
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

i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1; i_model = 1;
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
colors = {'r.','g.','b.','c.'};
MS_A = 12;
MS_B = 30; 
clf; hold on
axis square
low_lim  = -Inf;
high_lim = Inf;

if strcmp(details.metric,'BPS_divideCRM') || strcmp(details.metric,'BPS_divideUOP') || strcmp(details.metric, 'FracVar_10msec')
    low_lim  = 0;
    high_lim = 1;
end



if strcmp(details.metric, 'FracVar_10msec') || strcmp(details.metric, 'FracVar_20msec')
    low_lim  = 0;
    high_lim = 1;
end

if strcmp(details.metric,'VSPKD_50msec')
    low_lim = 0;
    high_lim = 2;
end
if strcmp(details.metric,'VSPKD_50msec_subtract')
    low_lim = 0;
    high_lim = 1;
end


if strcmp(details.metric,'VSPKD_50msec_divide')
    low_lim  = 1;
    high_lim = 2;
end



xlim([low_lim, high_lim]);
ylim([low_lim, high_lim]);
plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');
for i_exp = details.exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    colorstring = colors{i_exp};

    for i_celltype = details.celltypes
        model1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
        model2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;

        model1(find(model1<=low_lim)) = low_lim;
        model2(find(model2<=low_lim)) = low_lim;

        model1(find(model1>=high_lim)) = high_lim;
        model2(find(model2>=high_lim)) = high_lim;
        
        plot(model1, model2, colorstring, 'markersize', MS_B);
        if i_celltype == 1
            plot(model1, model2, 'w.', 'markersize', MS_A);
        elseif i_celltype == 2
            plot(model1, model2, 'k.', 'markersize', MS_A);
        end            
    end
end

cd(savedir)
eval(sprintf('print -dpdf plot_%s.pdf',savename))
cd(currentdir)    



%% PLOT THE 2 by 2, WN/NSEM and Fit vs Train
%{
if plot_diagnostics
    eval(sprintf('load %s/%s_model_comparsion_matfile_%s.mat model_comparison GLMType', savedir, comparison_name,cellselection_type))
    colors = {'r','g','b','c'};
    MS_A = 16;
    MS_B = 20;

    clf
    %{
    subplot(5,1,1)
    set(gca, 'fontsize', 10); axis off; c = 0;
    c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Type: %s', GLMType.fitname),'interpreter','none');
    c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none');
    subplot(5,1,[2 3 4 5]); hold on; set(gca,'fontsize',10)
    %}
    for i_exp = exps
        clear raster_scores
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
       % eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        colorstring = colors{i_exp};
        
        subplot(2,2,1); hold on;
        set(gca, 'fontsize', 10);
        title('WN: FIT: Percent Change in Objective Value');
        xlabel(sprintf('Raw Obj Val: %s', fit1),'interpreter','none');
        ylabel(sprintf('Pct Improve: %s', fit2),'interpreter','none');
        
        subplot(2,2,2); hold on
        set(gca, 'fontsize', 10);
        title('NSEM: FIT: Percent Change in Objective Value');
        xlabel(sprintf('Raw Obj Val: %s', fit1),'interpreter','none');
        ylabel(sprintf('Pct Improve: %s', fit2),'interpreter','none');
        
        
        subplot(2,2,3); hold on;
        set(gca, 'fontsize', 10);
        title('WN: CROSS VALIDATED: Fraction Bits');
        unity_line = linspace(0, 1,100);
        plot(unity_line,unity_line,'k')
        xlim([0,1]); xlabel(fit1,'interpreter','none')
        ylim([0,1]); ylabel(fit2,'interpreter','none')
        
        subplot(2,2,4); hold on
        set(gca, 'fontsize', 10);
        title('NSEM: CROSS VALIDATED: Fraction Bits');
        unity_line = linspace(0, 1,100);
        plot(unity_line,unity_line,'k')
        xlim([0,1]); xlabel(fit1,'interpreter','none')
        ylim([0,1]); ylabel(fit2,'interpreter','none')
        
        
        for i_celltype = celltypes
            xval_WN   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_bits;
            xval_NSEM = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_bits;
            xval_NSEM(find(xval_NSEM<0)) = 0;
            
            fit_WN   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_scores.WN.obj_val;
            fit_NSEM = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_scores.NSEM.obj_val;
                        
            if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
            if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
            plotstring_center = sprintf('%s.',colorstring);
            
            fitchange = (fit_WN(:,2)-fit_WN(:,1))./(fit_WN(:,1));
            %fitchange(find(fitchange<=0)) = 0;
         
            subplot(2,2,1); hold on;  
            plot(fit_WN(:,1), fitchange, plotstring_center,  'markersize', MS_B);
            plot(fit_WN(:,1), fitchange, plotstring_surround,  'markersize', MS_A);
            
            fitchange = (fit_NSEM(:,2)-fit_NSEM(:,1))./(fit_NSEM(:,1));
            %fitchange(find(fitchange<=0)) = 0;
            subplot(2,2,2); hold on;  
            plot(fit_NSEM(:,1), fitchange, plotstring_center,  'markersize', MS_B);
            plot(fit_NSEM(:,1), fitchange, plotstring_surround,  'markersize', MS_A);

            subplot(2,2,3); hold on;  
            plot(xval_WN(:,1), xval_WN(:,2), plotstring_center,  'markersize', MS_B);
            plot(xval_WN(:,1), xval_WN(:,2), plotstring_surround,  'markersize', MS_A);
            
            subplot(2,2,4); hold on
            plot(xval_NSEM(:,1), xval_NSEM(:,2), plotstring_center,  'markersize', MS_B);
            plot(xval_NSEM(:,1), xval_NSEM(:,2), plotstring_surround,  'markersize', MS_A);
        end

        
    end
    orient landscape
    eval(sprintf('print -dpdf %s/%s_WN_NSEM_FIT_XVAL_%s.pdf', savedir, comparison_name, cellselection_type))
end
%}


%% OLD CRAPPY PLOT TEMPLATES
 %{   
    figure
    colors = {'r','g','b','c'};
    MS_A = 16;
    MS_B = 20;
    minW = Inf;
    maxW = -Inf;
    minN = Inf;
    maxN = -Inf;
    %
    subplot(5,1,1)
    set(gca, 'fontsize', 10); axis off; c = 0;
    c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Type: %s', GLMType.fitname),'interpreter','none');
    c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none');
    subplot(5,1,[2 3 4 5]); hold on; set(gca,'fontsize',10)
    %
    for i_exp = exps
        clear raster_scores
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
       % eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        colorstring = colors{i_exp};

        for i_celltype = celltypes
            oldfit_WN   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_WN(:,1);
            oldfit_NSEM = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_NSEM(:,1);
            
            newfit_WN   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_WN(:,2);
            newfit_NSEM = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.fit_NSEM(:,2);
            
                        %scores_NSEM(find(scores_NSEM<-1)) = -1;
            
            if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
            if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
            plotstring_center = sprintf('%s.',colorstring);

            subplot(2,1,1); hold on;  
            plot(oldfit_WN, (newfit_WN-oldfit_WN)./(oldfit_WN), plotstring_center,  'markersize', MS_B);
            plot(oldfit_WN, (newfit_WN-oldfit_WN)./(oldfit_WN), plotstring_surround,  'markersize', MS_A);
            
            subplot(2,1,2); hold on;  
            plot(oldfit_NSEM, (newfit_NSEM-oldfit_NSEM)./(oldfit_NSEM), plotstring_center,  'markersize', MS_B);
            plot(oldfit_NSEM, (newfit_NSEM-oldfit_NSEM)./(oldfit_NSEM), plotstring_surround,  'markersize', MS_A);
            
            %
            minW = min(minW, min(scores_WN(:)));
            maxW = max(maxW, max(scores_WN(:)));
            
            minN = min(minN, min(scores_NSEM(:)));
            maxN = max(maxN, max(scores_NSEM(:)));
            %
            
        end

        
        %

    %
    subplot(2,1,1); hold on
        unity_line = linspace(minW, maxW,100);
        plot(unity_line,unity_line,'k')
        xlim([minW,maxW]);
        ylim([minW,maxW])
        
        subplot(2,1,2); hold on
        unity_line = linspace(minN, maxN,100);
        plot(unity_line,unity_line,'k')
        xlim([minN,maxN]);
        ylim([minN,maxN]);
%
    
    
%orient landscape
%eval(sprintf('print -dpdf %s/ALLrast_vsGLM_BPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))




%
eval(sprintf('save %s/model_comparison_%s.mat model_comparison', outputdir, GLMType.fit_type));
%%  Plot change in xval score
    clf;
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    if ~strcmp(Type, 'ConeModelComparison')
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    else
        text(0, 1-0.1*c,sprintf('Fit by: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    end
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit 1 on x-axis: %s',fittedGLM1.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit 2 on y axis: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))

    subplot(3,1,[2 3])
    plot(linspace(0,1,100), linspace(0,1,100),'k' ); hold on
    set(gca, 'fontsize', 14);
    set(gca,'xtick',[0:.2:.8]); set(gca,'ytick',[0:.2:.8]); 
    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        fit1   = model_comparison.byexpnm{i_exp}.normedbps_fit1;
        fit2   = model_comparison.byexpnm{i_exp}.normedbps_fit2;  

        a = find(fit1<0); fit1(a) = 0; 
        b = find(  fit2<0);   fit2(b) = 0; 

        maxval= max( maxval , max(max(fit1),max(fit2)) );
        minval= min( minval , min(max(fit1),max(fit2)) );
        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot(fit1(ONP), fit2(ONP),'r.','markersize',MS ); end
        if i_exp == 2, plot(fit1(ONP), fit2(ONP),'g.','markersize',MS ); end
        if i_exp == 3, plot(fit1(ONP), fit2(ONP),'b.','markersize',MS ); end
        if i_exp == 4, plot(fit1(ONP), fit2(ONP),'c.','markersize',MS ); end 

        if i_exp == 1, plot(fit1(OFFP), fit2(OFFP),'r*','markersize',MS ); end
        if i_exp == 2, plot(fit1(OFFP), fit2(OFFP),'g*','markersize',MS ); end
        if i_exp == 3, plot(fit1(OFFP), fit2(OFFP),'b*','markersize',MS ); end
        if i_exp == 4, plot(fit1(OFFP), fit2(OFFP),'c*','markersize',MS ); end 
    end


    xlim([0,1.1*maxval]);
    ylim([0,1.1*maxval]); hold on;
    title('Normalized Bits per Spike')
    xlabel(model_comparison.xlabel, 'interpreter','none')
    ylabel(model_comparison.ylabel,'interpreter','none')    
    eval(sprintf('print -dpdf %s/XValScores_%s.pdf', outputdir,GLMType.fit_type));


%
subplot(3,2,[4 6])
plot(linspace(0,1,100), linspace(0,1,100),'k' ); hold on
set(gca, 'fontsize', 14);
set(gca,'xtick',[0:.2:.8]); set(gca,'ytick',[0:.2:.8]); 
MS = 26;
maxval = 0;
minval = 0;
for i_exp = 1:4
    fit1   = model_comparison.byexpnm{i_exp}.normedbps_fit1;
    fit2   = model_comparison.byexpnm{i_exp}.normedbps_fit2;  
    
    a = find(fit1<=0);  
    b = find(  fit2<=0); 
    bad = union(a,b);
    
    param_vec = (fit2 - fit1) ./fit1 ;
    param_vec(bad) = 0;
    
    maxval = max(maxval, max(param_vec));
    minval = min(minval, min(param_vec));
    ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
    OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
    
    
    
    
    ONP  = setdiff(ONP, bad);
    OFFP = setdiff(OFFP,bad);
    if i_exp == 1, plot(1*ones(1,length(ONP)), param_vec(ONP), 'r.','markersize',MS ); end
    if i_exp == 2, plot(2*ones(1,length(ONP)), param_vec(ONP), 'g.','markersize',MS ); end
    if i_exp == 3, plot(3*ones(1,length(ONP)), param_vec(ONP), 'b.','markersize',MS ); end
    if i_exp == 4, plot(4*ones(1,length(ONP)), param_vec(ONP), 'c.','markersize',MS ); end 
        
    if i_exp == 1, plot(1*ones(1,length(OFFP)), param_vec(OFFP), 'r*','markersize',MS ); end
    if i_exp == 2, plot(2*ones(1,length(OFFP)), param_vec(OFFP), 'g*','markersize',MS ); end
    if i_exp == 3, plot(3*ones(1,length(OFFP)), param_vec(OFFP), 'b*','markersize',MS ); end
    if i_exp == 4, plot(4*ones(1,length(OFFP)), param_vec(OFFP), 'c*','markersize',MS ); end     
end


xlim([0 5]);
ylim([minval,1.1*maxval]); hold on;
title('Percent Improvement in Normalized Bits Per Spike')
xlabel(model_comparison.xlabel, 'interpreter','none')
ylabel(model_comparison.ylabel,'interpreter','none')  
%



%% show param clustering of 1-D double opt

if doubleopt  && (doubleopt_dim == 1)
    clf;

    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(Type , 'Stim_Nonlinearity')
        text(0, 1-0.1*c,sprintf('Black line or is the parameter values at which we there is no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black line is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(Type , 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
            param_meaning = 'Ratio of rescaling of values above the mean to values below the mean';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'raisepower_meanafter')
            param_meaning = 'Power to which we raised the input stimulus (orig on 0 to 1 scale)';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to HIGHER stimulus values';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to LOWER  stimulus values';
        end

        if strcmp(GLMType2.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            param_meaning = 'Power to which we raised the input stimulus about mean made odd ';
            above_statement = 'Values ABOVE black line indicate we give comparatively more weight to extreme highs and extreme lows of the stimulus values';
            below_statement = 'Values BELOW black line indicate we give comparatively more weight to middle values of the stimulus value';
        end
    elseif strcmp(Type , 'FilterOutput_Nonlinearity')
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'piece_linear_aboutmean')
            param_meaning = 'Ratio of rescaling of values of stimulus driven log CIF which are positive and negative';
            above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
        end
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            param_meaning = 'Ratio of rescaling of values of stimulus driven log CIF which are positive and negative';
            above_statement = 'Values ABOVE black line give added weight to extreme values (pos and neg) in lcif';
            below_statement = 'Values BELOW black line indicate gives added weight to values near 0';
        end    
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter (y-axis): %s', param_meaning))



    subplot(3,1,[2 3])
    set(gca, 'fontsize', 14); hold on


    xlim0 = 0; xlim1 = 5;
    if strcmp(Type, 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
            linear_point = 1;
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'raisepower_meanafter')
            linear_point = 1;
        end
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            linear_point = 1;
        end
    end

    if strcmp(Type, 'FilterOutput_Nonlinearity')
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'piece_linear_aboutmean')
            linear_point = 1;
        end
        if strcmp(GLMType2.postfilter_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
            linear_point = 1;
        end
    end


    LW = 2;    
    plot(linspace(xlim0,xlim1,100), linear_point*ones(1,100), 'k', 'linewidth', LW);
    xlim([xlim0,xlim1]);



    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.pt_nonlinearity_param(:)])
        maxval= max( maxval , max(param_vec) );
        minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot(1*ones(1,length(ONP)), param_vec(ONP), 'r.','markersize',MS ); end
        if i_exp == 2, plot(2*ones(1,length(ONP)), param_vec(ONP), 'g.','markersize',MS ); end
        if i_exp == 3, plot(3*ones(1,length(ONP)), param_vec(ONP), 'b.','markersize',MS ); end
        if i_exp == 4, plot(4*ones(1,length(ONP)), param_vec(ONP), 'c.','markersize',MS ); end 

        if i_exp == 1, plot(1*ones(1,length(OFFP)), param_vec(OFFP), 'r*','markersize',MS ); end
        if i_exp == 2, plot(2*ones(1,length(OFFP)), param_vec(OFFP), 'g*','markersize',MS ); end
        if i_exp == 3, plot(3*ones(1,length(OFFP)), param_vec(OFFP), 'b*','markersize',MS ); end
        if i_exp == 4, plot(4*ones(1,length(OFFP)), param_vec(OFFP), 'c*','markersize',MS ); end          
    end
    set(gca,'xtick',[1 2 3 4]);
    set(gca,'xticklabel', {'expA', 'expB', 'expC', 'expD'})
    if strcmp(Type, 'Stim_Nonlinearity')
        title('Optimized Parameter of Stimulus Non-linearity')
    end
    if strcmp(Type, 'FilterOutput_Nonlinearity')
        title('Optimized Parameter of Filter Output Non-linearity')
    end
    ylabel('Non-linearity Parameter')    
    orient landscape
    eval(sprintf('print -dpdf %s/OptimizedParamClustering_%s.pdf', outputdir,GLMType.fit_type));
end


%% show param clustering for 2-D double opt
if doubleopt && doubleopt_dim == 2
    clf
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(Type , 'Stim_Nonlinearity')
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value with no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(Type , 'Stim_Nonlinearity')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
            x_param_meaning = 'Ratio of rescaling of values above center to values below the mean';
            x_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            x_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            x_label           = 'Inc to Dec Ratio'
            
            y_param_meaning = 'Shift of Center away from additive mean';
            y_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            y_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            
            y_label          = 'New Center';
            
            no_modulation = [1,0];
            xlim0 = 0; xlim1 = 2;
            ylim0 = -.25; ylim1 = .25; 
        end  
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(x-axis): %s', x_param_meaning));
    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(y-axis): %s', y_param_meaning));
    
    subplot(3,1,[2 3]); xlim([xlim0,xlim1]); ylim([ylim0,ylim1]); hold on
    plot(no_modulation(1), no_modulation(2), 'k.', 'markersize',40); 
    plot(linspace(xlim0,xlim1,100), no_modulation(2)*ones(1,100),'k');
    plot(no_modulation(1)*ones(1,100),linspace(ylim0,ylim1,100), 'k');

    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.pt_nonlinearity_param(:)]);
        maxval= max( maxval , max(param_vec) );
        minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot( param_vec(ONP,1), param_vec(ONP,2), 'r.','markersize',MS ); end
        if i_exp == 2, plot( param_vec(ONP,1), param_vec(ONP,2), 'g.','markersize',MS ); end
        if i_exp == 3, plot( param_vec(ONP,1), param_vec(ONP,2), 'b.','markersize',MS ); end
        if i_exp == 4, plot( param_vec(ONP,1), param_vec(ONP,2), 'c.','markersize',MS ); end 

        if i_exp == 1, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'r*','markersize',MS ); end
        if i_exp == 2, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'g*','markersize',MS ); end
        if i_exp == 3, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'b*','markersize',MS ); end
        if i_exp == 4, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'c*','markersize',MS ); end          
    end
    if strcmp(Type, 'Stim_Nonlinearity')
        title('Optimized Parameter of Stimulus Non-linearity')
    end
    if strcmp(Type, 'FilterOutput_Nonlinearity')
        title('Optimized Parameter of Filter Output Non-linearity')
    end
    ylabel(y_label); xlabel(x_label)    
    orient landscape
    eval(sprintf('print -dpdf %s/OptimizedParamClustering_%s.pdf', outputdir,GLMType.fit_type));
end
%%
if doubleoptcomp 
    clf
    subplot(3,1,1);
    axis off
    set(gca, 'fontsize', 12)
    c = 1;
    text(0, 1-0.1*c,sprintf('%s', model_comparison.title),'interpreter','none' )
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit by: %s, Cone Model: %s',GLMType.fit_type, GLMType.cone_model),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Fit name: %s',fittedGLM2.GLMType.fitname),'interpreter','none')
    c = c+1;
    text(0, 1-0.1*c,sprintf('Dots are On-Parsols, Asterisks are Off-Parasols'))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Each Color is a different experiment'))
    c = c+1;
    if strcmp(modification,'inputNL')
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value with no stimulus modulation'))
    elseif strcmp(Type , 'FilterOutput_Nonlinearity') 
        text(0, 1-0.1*c,sprintf('Black DOT is the parameter value at which we there is no modulation of filter outputs'))
    end


    if strcmp(modification,'inputNL')
        if strcmp(GLMType2.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
            x_param_meaning = 'Ratio of rescaling of values above center to values below the mean';
            x_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            x_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            x_label           = 'Inc to Dec Ratio'
            
            y_param_meaning = 'Shift of Center away from additive mean';
            y_above_statement = 'Values ABOVE black line indicate comporatively more weight given to stim values ABOVE the original mean';
            y_below_statement = 'Values BELOW black line indicate comporatively more weight given to stim values BELOW the original mean';
            
            y_label          = 'New Center';
            
            no_modulation = [0,0];
            xlim0 = -1; xlim1 = 1;
            ylim0 = -.3; ylim1 = .3; 
        end  
    end

    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(x-axis): %s', x_param_meaning));
    c = c+1;
    text(0, 1-0.1*c,sprintf('Nonlinearity Parameter(y-axis): %s', y_param_meaning));
    
    subplot(3,1,[2 3]); xlim([xlim0,xlim1]); ylim([ylim0,ylim1]); hold on
    %plot(no_modulation(1), no_modulation(2), 'k.', 'markersize',40); 
    plot(linspace(xlim0,xlim1,100), no_modulation(2)*ones(1,100),'k');
    plot(no_modulation(1)*ones(1,100),linspace(ylim0,ylim1,100), 'k');

    MS = 26;
    maxval = 0;
    minval = 1;
    for i_exp = 1:4
        param_vec = cell2mat([model_comparison.byexpnm{i_exp}.param_search2(:)]) - cell2mat([model_comparison.byexpnm{i_exp}.param_search1(:)])
      %  maxval= max( maxval , max(param_vec) );
       % minval= min( minval , min(param_vec) );

        ONP  = find(model_comparison.byexpnm{i_exp}.ONP);
        OFFP = find(model_comparison.byexpnm{i_exp}.OFFP);
        if i_exp == 1, plot( param_vec(ONP,1), param_vec(ONP,2), 'r.','markersize',MS ); end
        if i_exp == 2, plot( param_vec(ONP,1), param_vec(ONP,2), 'g.','markersize',MS ); end
        if i_exp == 3, plot( param_vec(ONP,1), param_vec(ONP,2), 'b.','markersize',MS ); end
        if i_exp == 4, plot( param_vec(ONP,1), param_vec(ONP,2), 'c.','markersize',MS ); end 

        if i_exp == 1, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'r*','markersize',MS ); end
        if i_exp == 2, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'g*','markersize',MS ); end
        if i_exp == 3, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'b*','markersize',MS ); end
        if i_exp == 4, plot( param_vec(OFFP,1), param_vec(OFFP,2), 'c*','markersize',MS ); end          
    end
    title('Manual Params - FMINCON Params')
    ylabel(y_label); xlabel(x_label)    
    orient landscape
    eval(sprintf('print -dpdf %s/SearchDifferences_%s.pdf', outputdir,GLMType.fit_type));
end
end


%}


% Other Calls
%{
comparison_name = 'InputNL_Hinge_spatailSTA_linearCones';
glmsettings_x{1}.type = 'cone_model';
glmsettings_x{1}.name = 'rieke_linear';
glmsettings_y{1}.type = 'cone_model';
glmsettings_y{1}.name = 'rieke_linear';
glmsettings_y{2}.type = 'input_pt_nonlinearity';
glmsettings_y{2}.name = 'piece_linear_aboutmean';
fit1 = 'spatialSTA_linearCones';
fit2 = 'Hinge';


comparison_name = 'Convergence_Tol2_vs_Tol5_spatialSTA';
glmsettings_x = {};
glmsettings_y{1}.type = 'specialchange'
glmsettings_y{1}.name = 'TolFun_2'
fit1 = 'spatialSTA Tol5';
fit2 = 'spatialSTA Tol2';

comparison_name = 'Convergence_Tol3_vs_Tol5_spatialSTA';
glmsettings_x = {};
glmsettings_y{1}.type = 'specialchange'
glmsettings_y{1}.name = 'TolFun_3'
fit1 = 'spatialSTA Tol5';
fit2 = 'spatialSTA Tol3';


comparison_name = 'Rank1_vs_spatialSTA';
glmsettings_x = {};
glmsettings_y{1}.type = 'filter_mode';
glmsettings_y{1}.name = 'rk1';
fit1 = 'spatialSTA';
fit2 = 'Rank1';



comparison_name = 'Rank2_vs_Rank1';
glmsettings_x{1}.type = 'filter_mode';
glmsettings_x{1}.name = 'rk1';
glmsettings_y{1}.type = 'filter_mode';
glmsettings_y{1}.name = 'rk2';
fit1 = 'Rank1';
fit2 = 'Rank2';

comparison_name = 'Rank2_vs_spatialSTA';
glmsettings_x = [];
glmsettings_y{1}.type = 'filter_mode';
glmsettings_y{1}.name = 'rk2';
fit1 = 'spatialSTA';
fit2 = 'Rank2';

comparison_name = 'Rank2_vs_spatialSTA_LinearCones';
glmsettings_x{1}.type = 'cone_model';
glmsettings_x{1}.name = 'rieke_linear';
glmsettings_y{1}.type = 'cone_model';
glmsettings_y{1}.name = 'rieke_linear';
glmsettings_y{2}.type = 'filter_mode';
glmsettings_y{2}.name = 'rk2';
fit1 = 'spatialSTA';
fit2 = 'Rank2';

comparison_name = 'hinge_LinearCones';
glmsettings_x{1}.type = 'cone_model';
glmsettings_x{1}.name = 'rieke_linear';
glmsettings_y{1}.type = 'cone_model';
glmsettings_y{1}.name = 'rieke_linear';
glmsettings_y{2}.type = 'input_pt_nonlinearity';
glmsettings_y{2}.name = 'piece_linear_aboutmean';

comparison_name = 'InputNL_4piecelin';
glmsettings_x{1}.type = 'cone_model';
glmsettings_x{1}.name = 'rieke_linear';
glmsettings_y{1}.type = 'cone_model';
glmsettings_y{1}.name = 'rieke_linear';
glmsettings_y{2}.type= 'input_pt_nonlinearity';
glmsettings_y{2}.name= 'piecelinear_fourpiece_eightlevels';
fit1 = 'spatialSTA';
fit2 = '4piecelin';
%}

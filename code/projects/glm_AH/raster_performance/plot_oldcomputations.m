% AKHEITMAN 2015-03-09 START 
% NEW PLOT
% STRAIGHT FOR GLM_PERFORMANCE
clear; close all; clc
conv_type   = 'glmconv_4pct';
%conv_type   = 'glmconv_1pct';
%{
%m{etric_type = 'vspkd_logdivide_64bins';
%metric_type = 'vspkd_divide_64bins';
%metric_type = 'vspkd_divide_16bins';
metric_type = 'vspkd_subtract_64bins';
%metric_type = 'vspkd_subtract_16bins';
eval(sprintf('load Viktor_Spike_holdonto_2015-03-09.mat'))
%}

%
metric_type = 'fracvar_13bins';
%metric_type = 'fracvar_65bins';
if strcmp(metric_type, 'fracvar_13bins') || strcmp(metric_type, 'fracvar_65bins')
    min_min = 0;
end
eval(sprintf('load fracvar_avgsignals_2015-03-09.mat'))
%}


% RESET
if strcmp(conv_type,'glmconv_4pct')
    conv_column = 2;   
end
if strcmp(conv_type, 'glmconv_1pct')
    conv_column = 4;
end

% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE
celltypes = [1 2];
exps      = [1 2 3]; 

%{
debug = false;
compute = false; 
plotrasterscores = true;  
comparetoGLM = true;
%}

% SET DIRECTORIES / GLM WITH TIME KERNEL
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection));
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType);
savedir = sprintf('%s/PrototypePlots/WN_vs_NSEM_withGLMCONV/%s', BD.NSEM_home,GLMType.fitname);
if ~exist(savedir,'dir'), mkdir(savedir); end
% crm_type = 'base_crm_findPS';
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%% INDIVIDUAL PREPS

colors = {'r','g','b','c'};
MS_A = 12;
MS_B = 16;

for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    colorstring = colors{i_exp};

    for i_celltype = celltypes
        clf;
        if i_celltype == 1, celltype = 'ONP';  cellgroup0 = allcells{i_exp}.ONP; conv_index = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup0 = allcells{i_exp}.OFFP;conv_index = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column)); end
        
        cellgroup = cellgroup0(conv_index);
        % NEW COMPUTATION
        %{
        scores_WN   = scores.celltype{i_celltype}.scores_WN;
        scores_NSEM = scores.celltype{i_celltype}.scores_NSEM;
        %}
        
        % OLD STRUCTURE
        if strcmp(metric_type,'vspkd_divide_64bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(5,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(5,conv_index);
        end
        if strcmp(metric_type,'vspkd_logdivide_64bins')
            scores_WN   = log(glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(5,conv_index));
            scores_NSEM = log(glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(5,conv_index));
        end
        if strcmp(metric_type,'vspkd_divide_16bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(3,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(3,conv_index);
        end
        if strcmp(metric_type,'vspkd_subtract_64bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract(5,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract(5,conv_index);
        end
        if strcmp(metric_type,'vspkd_subtract_16bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract(3,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract(3,conv_index);
        end
        
        if strcmp(metric_type,'fracvar_13bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw(7,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw(7,conv_index);
        end
        if strcmp(metric_type,'fracvar_65bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw(14,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw(14,conv_index);
        end
        
        if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
        if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
        plotstring_center = sprintf('%s.',colorstring);
        
        subplot(5,1,1)
        set(gca, 'fontsize', 10); axis off
        c = 0;
        c=c+1; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s, Metric: %s',exp_nm, celltype,metric_type),'interpreter','none');
        c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
        %c=c+1; text(-.1, 1-0.1*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
        %c=c+1; text(-.1, 1-0.1*c,sprintf('BPS: Bits Per Spike - Improvement from Flat Rate Model, more positive is better'));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
    
        for i_metric = 1
            vals_WN   = scores_WN;
            vals_NSEM = scores_NSEM; 
            
            % Throw out Garbage
            bad_ind = union(find(abs(vals_WN) == Inf),find(abs(vals_NSEM) == Inf));
            vals_WN(bad_ind) = [];
            vals_NSEM(bad_ind) = [];
            
            
            if exist('min_min', 'var')
                vals_WN(find(vals_WN<=min_min)) = min_min;
                vals_NSEM(find(vals_NSEM<=min_min)) = min_min;
            end
            
            
            maxval = max(max(vals_WN),max(vals_NSEM));
            minval = min(min(vals_WN),min(vals_NSEM));
            
            
            
            % Throw out 
            subplot(5,1,[2,3,4,5]); hold on; set(gca,'fontsize',10)
            maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',metric_type),'interpreter','none')

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
            xlabel('White Noise');
            ylabel('Natural Scenes');
        end
        
        orient tall
        orient landscape
        eval(sprintf('print -dpdf %s/%s_%s_%s_%s.pdf', savedir, metric_type,conv_type, exp_nm, celltype))
    end
end


clf
maxval = -Inf;
minval = Inf;
subplot(5,1,1)
set(gca, 'fontsize', 10); axis off
c = 0;
c=c+1; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s, Metric: %s',exp_nm, celltype,metric_type),'interpreter','none');
c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
        %c=c+1; text(-.1, 1-0.1*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
        %c=c+1; text(-.1, 1-0.1*c,sprintf('BPS: Bits Per Spike - Improvement from Flat Rate Model, more positive is better'));
c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
for i_exp = exps
    subplot(5,1,[2,3,4,5]); hold on; set(gca,'fontsize',10)
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    colorstring = colors{i_exp};

    for i_celltype = celltypes
        if i_celltype == 1, celltype = 'ONP';  cellgroup0 = allcells{i_exp}.ONP; conv_index = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup0 = allcells{i_exp}.OFFP;conv_index = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column)); end
        
        cellgroup = cellgroup0(conv_index);
        % NEW COMPUTATION
        %{
        scores_WN   = scores.celltype{i_celltype}.scores_WN;
        scores_NSEM = scores.celltype{i_celltype}.scores_NSEM;
        %}
        
        % OLD STRUCTURE
        if strcmp(metric_type,'vspkd_divide_64bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(5,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(5,conv_index);
        end
         if strcmp(metric_type,'vspkd_logdivide_64bins')
            scores_WN   = log(glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(5,conv_index));
            scores_NSEM = log(glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(5,conv_index));
        end
        if strcmp(metric_type,'vspkd_divide_16bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(3,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(3,conv_index);
        end
        if strcmp(metric_type,'vspkd_subtract_64bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract(5,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract(5,conv_index);
        end
        if strcmp(metric_type,'vspkd_subtract_16bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract(3,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract(3,conv_index);
        end
        if strcmp(metric_type,'fracvar_13bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw(7,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw(7,conv_index);
        end
        if strcmp(metric_type,'fracvar_65bins')
            scores_WN   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw(14,conv_index);
            scores_NSEM = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw(14,conv_index);
        end
        
        if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
        if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
        plotstring_center = sprintf('%s.',colorstring);
        
        
    
        for i_metric = 1
            vals_WN   = scores_WN;
            vals_NSEM = scores_NSEM; 
            
            % Throw out Garbage
            bad_ind = union(find(abs(vals_WN) == Inf),find(abs(vals_NSEM) == Inf));
            vals_WN(bad_ind) = [];
            vals_NSEM(bad_ind) = [];
            if exist('min_min', 'var')
                vals_WN(find(vals_WN<=min_min)) = min_min;
                vals_NSEM(find(vals_NSEM<=min_min)) = min_min;
            end

            newmax = max(max(vals_WN),max(vals_NSEM));
            newmin = min(min(vals_WN),min(vals_NSEM));
            
            maxval = max(newmax,maxval);
            minval = min(newmin,minval);

            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',metric_type),'interpreter','none')

            if exist('label_points','var') && label_points
                for i_cell = 1:length(cids)
                    lab = sprintf('%d', cids(i_cell));
                    text(vals(1,i_cell),vals(2,i_cell),lab);
                end
            end
            
        end
        
        
    end
end
unity_line = linspace(minval, maxval,100);
plot(unity_line,unity_line,'k')
xlim([minval,maxval]);
ylim([minval,maxval]);
xlabel('White Noise');
ylabel('Natural Scenes');

orient tall
orient landscape
eval(sprintf('print -dpdf %s/%s_%s_ALLTOGETHER.pdf', savedir, metric_type,conv_type))
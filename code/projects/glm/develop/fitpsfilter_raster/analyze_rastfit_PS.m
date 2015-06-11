% RESET
% 2015-05-31
clear; close all; clc

% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE
celltypes = [1 2]
exps = [1 2 3] 
analyze_PSgain = true;  
analyze_PSplot = true;

comparetoGLM_withCONV = false;

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
savedir = sprintf('%s/Raster_Metrics', BD.Cell_Selection);
if ~exist(savedir,'dir'), mkdir(savedir); end

% CRM = CONDITIONALIZED RATE MODEL
%crm_type = 'base_crm_importPS';
crm_type = 'base_crm_findPS_nosignal';
%crm_type = 'base_crm_findPS_unitlikebasis';
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
%% PLOTTING RAW RASTER SCORES
if analyze_PSgain
% Plot WN-NSEM comparison for these base predictability metrics
colors = {'r','g','b','c'};
exp_strings = {'expA','expB','expC'};
clf;
subplot(4,1,1);
set(gca, 'fontsize', 8); axis off
c = 0; text(-.1, 1-0.15*c,sprintf('Analysing PS filter fit from Raster'), 'interpreter','none');
c=c+1; text(-.1, 1-0.15*c,sprintf('Plotting Geometric Mean of PS Gain')) ;
c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;

MS_OUT = 30;
MS_IN  = 10;

low_cut  = .5;
high_cut = 1;
for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
    colorstring = colors{i_exp};
    exp_string =  exp_strings{i_exp}; 
    for i_celltype = celltypes
        
        if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end
      
        ps_WN   = raster_scores.celltype{i_celltype}.scores_WN.ps_filter;
        ps_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM.ps_filter;
        
        plotstring_surround = sprintf('%s.', colorstring);
        if i_celltype == 1, plotstring_center = 'w.'; titlestring = 'ONP'; end
        if i_celltype == 2, plotstring_center = 'k.'; titlestring = 'OFFP';end
        
        vals_WN = NaN(1, length(ps_WN));
        for i_cell = 1:length(ps_WN)
            vals_WN(i_cell) = geomean(exp(ps_WN{i_cell}));
        end
        vals_NSEM = NaN(1, length(ps_NSEM));
        for i_cell = 1:length(ps_NSEM)
            vals_NSEM(i_cell) = geomean(exp(ps_NSEM{i_cell}));
        end
        subplot(4,2,(i_exp*2+i_celltype)); hold on; set(gca,'fontsize',8)
        plot(linspace(low_cut,high_cut,100),linspace(low_cut,high_cut,100),'k')
        
        vals_WN(find(vals_WN<=low_cut))      = low_cut;
        vals_NSEM(find(vals_NSEM<=low_cut))  = low_cut;
        vals_WN(find(vals_WN>=high_cut))     = high_cut;
        vals_NSEM(find(vals_NSEM>=high_cut)) = high_cut;
        
        plot(vals_WN, vals_NSEM, plotstring_surround,  'markersize', MS_OUT);
        plot(vals_WN, vals_NSEM, plotstring_center,'markersize', MS_IN);
        
        xlim([low_cut high_cut]);
        ylim([low_cut high_cut]);
        
        set(gca,'xtick', [.5 .75 1]);
        set(gca,'ytick', [.5 .75 1]);
        
        xlabel('White Noise');
        ylabel('Natural Scenes');
        title(sprintf('Mean Gain of Raster Fitted PS Filter: %s %s', exp_string, titlestring))

        orient tall
        eval(sprintf('print -dpdf %s/rast_fittedPS_meanVAL_%s.pdf', savedir, crm_type))
    end
end


end


%%
if analyze_PSplot

for i_fit = 1:2
    colors = {'r','g','b','c'};
    exp_strings = {'expA','expB','expC'};
    clf;subplot(4,1,1)
    set(gca, 'fontsize', 10); axis off; 
    c = 0; text(-.1, 1-0.15*c,sprintf('Overlay of all Raster Fit Post-Spike Filters, no cuts'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none');
    y_low_cut = 0;
    y_high_cut = 2;
    for i_exp = exps
        clear raster_scores
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        colorstring = colors{i_exp};
        exp_string =  exp_strings{i_exp};
        for i_celltype = celltypes
            if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;   end
            if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP;  end
            
            if i_fit == 1
                ps_filters   = raster_scores.celltype{i_celltype}.scores_WN.ps_filter;
                fit_type = 'WN';
            elseif i_fit == 2
                ps_filters = raster_scores.celltype{i_celltype}.scores_NSEM.ps_filter;
                fit_type = 'NSEM';
            end

            subplot(4,2,(i_exp*2+i_celltype)); hold on; set(gca,'fontsize',8)
            xlabel('post spike time [MSEC]');
            ylabel('Gain');
            title(sprintf('%s Raster Fitted PS Filters: %s %s', fit_type, exp_string, celltype)); 

            xlim([0 100]);
            ylim([y_low_cut y_high_cut]);
            set(gca,'ytick', [0 .5 1 1.5 2]);
            set(gca,'xtick', [0 25 50 75 100]);
            plot(linspace(0,100,120), ones(1,120), 'k', 'linewidth', 2);
            for i_cell =1:length(ps_filters)
                plot(linspace(0,100,120),exp(ps_filters{i_cell}), colorstring,'linewidth',.5);
            end
            
            
            
        end

    end
    orient tall
    eval(sprintf('print -dpdf %s/rast_fittedPS_plotOVERLAY_%s_%s.pdf', savedir, fit_type,crm_type))
end
end

%%
%{
conv_type   = 'glmconv_4pct';
if strcmp(conv_type,'glmconv_4pct')
    conv_column = 2;   
end


if comparetoGLM_withCONV
    colors = {'r','g','b','c'};
    MS_A = 12;
    MS_B = 16;
    subplot(5,1,1)
    set(gca, 'fontsize', 10); axis off; c = 0;
    c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
    c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Type: %s', GLMType.fitname),'interpreter','none');
    c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none');
    subplot(5,1,[2 3 4 5]); hold on; set(gca,'fontsize',10)
    for i_exp = exps
        clear raster_scores
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        colorstring = colors{i_exp};

        for i_celltype = celltypes
            
            
            if i_celltype == 1, celltype = 'ONP';  cellgroup0 = allcells{i_exp}.ONP;  conv_index = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column)); end
            if i_celltype == 2, celltype = 'OFFP'; cellgroup0 = allcells{i_exp}.OFFP; conv_index = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));end
            cellgroup = cellgroup0(conv_index);
            
            scores_WN   = raster_scores.celltype{i_celltype}.scores_WN;
            scores_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM;
            
            
            

            if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
            if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
            plotstring_center = sprintf('%s.',colorstring);

            glm_WN    = scores_WN.glm_bps(conv_index); 
            glm_NSEM  = scores_NSEM.glm_bps(conv_index);   
            norm_WN   = scores_WN.crm_bps(conv_index); 
            norm_NSEM = scores_NSEM.crm_bps(conv_index);
            vals_WN = glm_WN ./ norm_WN;
            vals_NSEM = glm_NSEM ./norm_NSEM;
            
            bad_ind = union(find(abs(vals_WN) == Inf),find(abs(vals_NSEM) == Inf));
            vals_WN(bad_ind) = [];
            vals_NSEM(bad_ind) = [];
            
            vals_WN(find(vals_WN<=0))= 0;
            vals_NSEM(find(vals_NSEM<=0)) = 0;
            
            
                maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
                minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
                plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
                plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
                title(sprintf('%s','BPS GLM / BPS Conditioned Rate Model'))

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
    xlabel('White Noise');
    ylabel('Natural Scenes');
end
orient landscape
eval(sprintf('print -dpdf %s/ALLrast_vsGLM_BPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))
%}




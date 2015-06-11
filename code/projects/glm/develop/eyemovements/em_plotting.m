clear ; close all; clc;
cellselection_type = 'glmconv_4pct';
%cellselection_type = 'shortlist';
load_data = true;
plot_diagnostics = false;
plot_deltaWNNSEM = true;

comparison_name = 'overfitted_fem_saccade';
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear';
changes_cell{2}.type= 'input_pt_nonlinearity';
changes_cell{2}.name= 'piecelinear_fourpiece_eightlevels';
fit1 = 'GLM';
fit2 = 'withEM'

celltypes = [1 2];
exps = [1 2 3];
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
savedir = sprintf('%s/PrototypePlots/ModelComparisons/%s',BD.NSEM_home,comparison_name);
if ~exist(savedir,'dir'), mkdir(savedir); end
GLMType = GLM_settings('default',changes_cell);
GLMType.fitname    = GLM_fitname(GLMType);

%%% LOADING 
model_comparison.byexpnm         = allcells;
model_comparison.comparison_name = comparison_name;
model_comparison.fitname         = GLMType.fitname;
%model_comparison.normalizer      = normalizer;
model_comparison.fullGLMType     = GLMType;
model_comparison.notes.n1        = 'Use Normalization from ';
model_comparison.notes.n2        = 'UOP: unconditioned optimal performance';
model_comparison.notes.n3        = 'CRM: conditioned rate model';
model_comparison.notes.n4        = 'BPS: Bits Per Spike';
model_comparison.notes.n5        = 'LPPS: logarithmic probability per second';
model_comparison.notes.n6        = 'GLM: Generalized Linear Model';
model_comparison.timestamp       = datestr(clock);
model_comparison.code_name       = mfilename('fullpath');
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1; i_model = 1;

%% LOAD DATA FROM FITTED GLM
if load_data
    
% First Load BPS
for i_exp = exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    if i_exp < 4
        normalizer = 'base_crm_findPS';
    elseif i_exp == 4
        normalizer = 'base_crm_importPS';
    end   
    eval(sprintf('load %s/Raster_Metrics/%s_%s.mat raster_scores', BD.Cell_Selection, normalizer, exp_nm));
    raster_bps  = raster_scores; clear raster_scores

    

    secondDir.exp_nm    = exp_nm;
    secondDir.stim_type = 'WN';
    secondDir.map_type  = 'mapPRJ';
    secondDir.fitname   = GLMType.fitname;
    basedir.WN = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
    crossval_dir.WN = sprintf('%s/CrossVal_RawMetrics',basedir.WN);
    clear secondDir
    
    BD = NSEM_BaseDirectories
    BD.GLM_output_raw = sprintf('%s/EM/saccade_NEGPOS_longfem_displacement_NEGPOS', BD.GLM_output_raw);
    secondDir.exp_nm    = exp_nm;
    secondDir.stim_type = 'NSEM';
    secondDir.map_type  = 'mapPRJ';
    secondDir.fitname   = GLMType.fitname;
    basedir.NSEM = NSEM_secondaryDirectories('savedir_GLMfit', secondDir,'',BD);
    crossval_dir.NSEM = sprintf('%s/CrossVal_RawMetrics',basedir.NSEM);
    clear secondDir 

    candidate_cells = [allcells{i_exp}.ONP , allcells{i_exp}.OFFP];
    if strcmp(cellselection_type, 'shortlist')
        [exp_nm,candidate_cells,expname]  = cell_list(i_exp,cellselection_type);
        candidate_cells = cell2mat(candidate_cells);
    elseif strcmp(cellselection_type,'glmconv_4pct')
        conv_column = 2; 
        conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
        conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
        candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
        clear conv_column conv_index_ON conv_index_OFF
    end
    
    %%
    for i_celltype = celltypes
        if i_celltype == 1, celltype = 'ONPar';  cellgroup_full = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup_full = allcells{i_exp}.OFFP; end
        
        cellgroup     = intersect(cellgroup_full,candidate_cells);
        if i_celltype == 1, model_comparison.byexpnm{i_exp}.ONP  = cellgroup; end
        if i_celltype == 2, model_comparison.byexpnm{i_exp}.OFFP = cellgroup; end
        
        xval_scores.WN.frac_bits     = NaN(length(cellgroup),2);       
        xval_scores.NSEM.frac_bits   = NaN(length(cellgroup),2);
        
        rastnorm.WN.bps           = raster_bps.celltype{i_celltype}.scores_WN.crm_bps;
        %%
        % Load Cell Specific Elements   Spikes and STA
        for i_cell = 1:length(cellgroup)
            clear crossval_rawmetrics cell_savename
            cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
            display(sprintf('LOADING %s %s', expname,cell_savename));
            norm_index = find(cellgroup_full == cid);
            
            eval(sprintf('load %s/%s.mat fittedGLM', basedir.WN, cell_savename))
            rawXVAL.WN.bps(1,1) = fittedGLM.xvalperformance.logprob_glm_bpspike;
            rawXVAL.WN.bps(1,2) = fittedGLM.xvalperformance.logprob_glm_bpspike;
            xval_scores.WN.frac_bits(i_cell,:)   = rawXVAL.WN.bps   / rastnorm.WN.bps(norm_index);
            
            eval(sprintf('load %s/%s.mat xvalperformance_NEW', basedir.NSEM, cell_savename))
          
            %%
            xval_scores.NSEM.frac_bits(i_cell,1) = xvalperformance_NEW.logprob_glm_OLD_bps_normCRM;
            xval_scores.NSEM.frac_bits(i_cell,2) = xvalperformance_NEW.logprob_glm_NEW_bps_normCRM;
            clear rawXVVAL 
        end
        
        model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores = xval_scores;
    end
end
eval(sprintf('save %s/%s_model_comparsion_matfile_%s.mat model_comparison GLMType', savedir, comparison_name,cellselection_type))
end

%% PLOT THE 2 by 2, WN/NSEM and Fit vs Train
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


%% PLOT CHANGE IN WN VS NSEM  (ONE SINGLE LARGE PLOT)
close all
if plot_deltaWNNSEM
eval(sprintf('load %s/%s_model_comparsion_matfile_%s.mat model_comparison GLMType', savedir, comparison_name,cellselection_type))
colors = {'r.','g.','b.','c.'};
MS_A = 12;
MS_B = 30; 
    
for i_metric = 1
    clear score_max score_min plot_max plot_min
    if i_metric == 1, met_type = 'frac_bits'; score_max = 1; score_min = 0; plot_max = 1; plot_min = 0; end
    if i_metric == 2, met_type = 'vspkd_ratio'; score_max = 5; score_min= 1;plot_max = 5; plot_min = 1; end
    if i_metric == 3, met_type = 'vspkd_diff'; score_max = 1.6; score_min= 0; plot_max = 1.6; plot_min = 0; end
    if i_metric == 4, met_type = 'frac_var_unexp';  score_max = 1; score_min = 0; plot_max = 1; plot_min = 0; end
    if i_metric == 5, met_type = 'frac_var_exp';  score_max = 1; score_min = 0; plot_max = 1; plot_min = 0; end
    figure(1); clf; hold on; axis square
    figure(2); clf; hold on; axis square
    figure(3); clf; hold on; axis square
    for i_exp = exps
        clear raster_scores
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        colorstring = colors{i_exp};
        
        for i_celltype = celltypes
            if i_metric == 1
                scores_WN     = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_bits;
                scores_NSEM   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_bits;
            elseif i_metric == 2
                scores_WN     = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.vspkd_ratio);
                scores_NSEM   = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.vspkd_ratio);
            elseif i_metric == 3
                scores_WN     = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.vspkd_diff);
                scores_NSEM   = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.vspkd_diff);
            elseif i_metric == 4
                scores_WN     = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_var;
                scores_NSEM   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_var;
            
            elseif i_metric == 5
                scores_WN     = 1-model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_var;
                scores_NSEM   = 1-model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_var;
            end
            scores_WN(find(scores_WN>score_max)) = score_max;
            scores_WN(find(scores_WN<score_min)) = score_min;
            scores_NSEM(find(scores_NSEM>score_max)) = score_max;
            scores_NSEM(find(scores_NSEM<score_min)) = score_min;
            
            figure(1);
            plot(scores_WN(:,1), scores_NSEM(:,1), colorstring, 'markersize', MS_B);
            for i_cell = 1:size(scores_WN,1)
                x_line = linspace(scores_WN(i_cell,1), scores_WN(i_cell,2),100);
                y_line = linspace(scores_NSEM(i_cell,1), scores_NSEM(i_cell,2),100);
                plot(x_line,y_line,'k');
            end
            plot(scores_WN(:,2), scores_NSEM(:,2), colorstring, 'markersize', MS_B);
            if i_celltype == 1
                plot(scores_WN(:,2), scores_NSEM(:,2), 'w.', 'markersize', MS_A);
            elseif i_celltype == 2
                plot(scores_WN(:,2), scores_NSEM(:,2), 'k.', 'markersize', MS_A);
            end
            
            figure(2);
            plot(scores_WN(:,1), scores_NSEM(:,1), colorstring, 'markersize', MS_B);
            if i_celltype == 1
                plot(scores_WN(:,1), scores_NSEM(:,1), 'w.', 'markersize', MS_A);
            elseif i_celltype == 2
                plot(scores_WN(:,1), scores_NSEM(:,1), 'k.', 'markersize', MS_A);
            end
            
            figure(3);
            plot(scores_WN(:,2), scores_NSEM(:,2), colorstring, 'markersize', MS_B);
            if i_celltype == 1
                plot(scores_WN(:,2), scores_NSEM(:,2), 'w.', 'markersize', MS_A);
            elseif i_celltype == 2
                plot(scores_WN(:,2), scores_NSEM(:,2), 'k.', 'markersize', MS_A);
            end
            
        end
    end
    figure(1)
    unity_line = linspace(plot_min, plot_max,100);
    plot(unity_line,unity_line,'k')
    xlim([plot_min, plot_max]); ylim([plot_min, plot_max]);
    set(gca,'xtick',plot_min+(plot_max-plot_min)*[0:.25:1]); 
    set(gca,'ytick',plot_min+(plot_max-plot_min)*[0:.25:1]);
    
    figure(2)
    unity_line = linspace(plot_min, plot_max,100);
    plot(unity_line,unity_line,'k')
    xlim([plot_min, plot_max]); ylim([plot_min, plot_max]);
    set(gca,'xtick',plot_min+(plot_max-plot_min)*[0:.25:1]); 
    set(gca,'ytick',plot_min+(plot_max-plot_min)*[0:.25:1]);
    
    figure(3)
    unity_line = linspace(plot_min, plot_max,100);
    plot(unity_line,unity_line,'k')
    xlim([plot_min, plot_max]); ylim([plot_min, plot_max]);
    set(gca,'xtick',plot_min+(plot_max-plot_min)*[0:.25:1]); 
    set(gca,'ytick',plot_min+(plot_max-plot_min)*[0:.25:1]);
    
    figure(1);
    eval(sprintf('print -dpdf %s/%s_deltaWNNSEM_%s_%s.pdf', savedir, comparison_name,met_type, cellselection_type))
    figure(2);
    eval(sprintf('print -dpdf %s/%s_oldWNNSEM_%s_%s.pdf', savedir, comparison_name,met_type, cellselection_type))
    figure(3);
    eval(sprintf('print -dpdf %s/%s_newWNNSEM_%s_%s.pdf', savedir, comparison_name,met_type, cellselection_type))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    for i_plot = 1:2 
        clf; hold on; axis square
        for i_exp = exps
            clear raster_scores
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            colorstring = colors{i_exp};

            for i_celltype = celltypes
                if i_metric == 1
                scores_WN     = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_bits;
                scores_NSEM   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_bits;
                elseif i_metric == 2
                    scores_WN     = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.frac_var;
                    scores_NSEM   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.frac_var;
                elseif i_metric == 3
                    scores_WN     = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.WN.vspkd_ratio);
                    scores_NSEM   = (model_comparison.byexpnm{i_exp}.celltype{i_celltype}.xval_scores.NSEM.vspkd_ratio);
                end


                scores_WN(find(scores_WN>score_max)) = score_max;
                scores_WN(find(scores_WN<score_min)) = score_min;

                scores_NSEM(find(scores_NSEM>score_max)) = score_max;
                scores_NSEM(find(scores_NSEM<score_min)) = score_min;

                plot(scores_WN(:,i_plot), scores_NSEM(:,i_plot), colorstring, 'markersize', MS_B);
                if i_celltype == 1
                    plot(scores_WN(:,i_plot), scores_NSEM(:,i_plot), 'w.', 'markersize', MS_A);
                elseif i_celltype == 2
                    plot(scores_WN(:,i_plot), scores_NSEM(:,i_plot), 'k.', 'markersize', MS_A);
                end
            end
        end
        unity_line = linspace(plot_min, plot_max,100);
        plot(unity_line,unity_line,'k')
        xlim([plot_min, plot_max]); ylim([plot_min, plot_max]);
        if i_metric == 1 || i_metric ==2 
            set(gca,'xtick',[0:.25:1]); set(gca,'ytick',[0:.25:1]);
        end
        if i_plot == 1
            eval(sprintf('print -dpdf %s/%s_oldWNNSEM_%s_%s.pdf', savedir, comparison_name, met_type, cellselection_type))
        elseif i_plot == 2
            eval(sprintf('print -dpdf %s/%s_newWNNSEM_%s_%s.pdf', savedir, comparison_name,met_type, cellselection_type))
        end
    end
    %}
end

end


% Version 0 .. Hack.. NO distinction between fit and XVAL
% hack_comparemodels
% AKHEITMAN 2015-03-13
% 
% GOAL:
% Modernization of glm_AH/newtestcode/comparemodels.m
%
% CALLS:
% GLM_settings, GLM_fitname NSEM_BaseDirectories cell_list


clear ; close all; clc;
%cellselection_type = 'all';
cellselection_type = 'shortlist';
normalizer = 'base_crm_importPS';
celltypes = [1 2];
exps = [1 2 3 4];
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
savedir = sprintf('%s/PrototypePlots/Model_Comparisons', BD.Cell_Selection);
if ~exist(savedir,'dir'), mkdir(savedir); end

%
comparison_name = 'Rank2_vs_spatialSTA';
changes_cell_A = [];
changes_cell_B{1}.type = 'filter_mode';
changes_cell_B{1}.name = 'rk2';
%}

%{
comparison_name = 'Rank2_vs_spatialSTA_LinearCones';
changes_cell_A{1}.type = 'cone_model';
changes_cell_A{1}.name = 'rieke_linear';
changes_cell_B{1}.type = 'cone_model';
changes_cell_B{1}.name = 'rieke_linear';
changes_cell_B{2}.type = 'filter_mode';
changes_cell_B{2}.name = 'rk2';
%}

%{
comparison_name = 'hinge_LinearCones';
changes_cell_A{1}.type = 'cone_model';
changes_cell_A{1}.name = 'rieke_linear';
changes_cell_B{1}.type = 'cone_model';
changes_cell_B{1}.name = 'rieke_linear';
changes_cell_B{2}.type = 'input_pt_nonlinearity';
changes_cell_B{2}.name = 'piece_linear_aboutmean';
%}


GLMType{1} = GLM_settings('default',changes_cell_A);
GLMType{1}.fitname    = GLM_fitname(GLMType{1});
GLMType{2} = GLM_settings('default',changes_cell_B);
GLMType{2}.fitname    = GLM_fitname(GLMType{2});

%%% LOADING 
model_comparison.byexpnm         = allcells;
model_comparison.comparison_name = comparison_name;
model_comparison.fitname1        = GLMType{1}.fitname;
model_comparison.fitname2        = GLMType{2}.fitname;
model_comparison.normalizer      = normalizer;
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

%% Load up parameters and xval scores across cells and experiements
for i_exp = exps
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load %s/Raster_Metrics/%s_%s.mat raster_scores', BD.Cell_Selection, normalizer, exp_nm));
    
    
    for i_model = 1:2
        secondDir.exp_nm    = exp_nm;
        secondDir.stim_type = 'WN';
        secondDir.map_type  = 'mapPRJ';
        secondDir.fitname   = GLMType{i_model}.fitname;
        GLM_fitdir{i_model}.WN = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        clear secondDir
        secondDir.exp_nm    = exp_nm;
        secondDir.stim_type = 'NSEM';
        secondDir.map_type  = 'mapPRJ';
        secondDir.fitname   = GLMType{i_model}.fitname;
        GLM_fitdir{i_model}.NSEM = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        clear secondDir
    end
    
    candidate_cells = [allcells{i_exp}.ONP , allcells{i_exp}.OFFP];
    if strcmp(cellselection_type, 'shortlist')
        [exp_nm,candidate_cells,expname]  = cell_list(i_exp,cellselection_type);
        candidate_cells = cell2mat(candidate_cells);
    end
    
    %%
    for i_celltype = celltypes
        if i_celltype == 1, celltype = 'ONPar';  cellgroup_full = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup_full = allcells{i_exp}.OFFP; end
        
        cellgroup     = intersect(cellgroup_full,candidate_cells);
        if i_celltype == 1, model_comparison.byexpnm{i_exp}.ONP  = cellgroup; end
        if i_celltype == 2, model_comparison.byexpnm{i_exp}.OFFP = cellgroup; end
        
        scores_WN     = NaN(length(cellgroup),2);
        scores_NSEM   = NaN(length(cellgroup),2);
        norm_WN       = raster_scores.celltype{i_celltype}.scores_WN.crm_bps;
        norm_NSEM     = raster_scores.celltype{i_celltype}.scores_NSEM.crm_bps;
        
        % Load Cell Specific Elements   Spikes and STA
        for i_cell = 1:length(cellgroup)
            clear fittedGLM cell_savename
            cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
            rawWN   = NaN(1,2);
            rawNSEM = NaN(1,2);
            for i_model = 1:2
                eval(sprintf('load %s/%s.mat fittedGLM', GLM_fitdir{i_model}.WN, cell_savename))
                rawWN(i_model) = fittedGLM.xvalperformance.logprob_glm_bpspike;
                clear fittedGLM
                eval(sprintf('load %s/%s.mat fittedGLM', GLM_fitdir{i_model}.NSEM, cell_savename))
                rawNSEM(i_model) = fittedGLM.xvalperformance.logprob_glm_bpspike;
                clear fittedGLM
            end
            scores_WN(i_cell,:)   = rawWN;
            scores_NSEM(i_cell,:) = rawNSEM;
        end
        
        model_comparison.byexpnm{i_exp}.celltype{i_celltype}.scores_WN = scores_WN;
        model_comparison.byexpnm{i_exp}.celltype{i_celltype}.scores_NSEM = scores_NSEM;
    end
end

%%
figure
    colors = {'r','g','b','c'};
    MS_A = 16;
    MS_B = 20;
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

        for i_celltype = celltypes
            scores_WN   = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.scores_WN;
            scores_NSEM = model_comparison.byexpnm{i_exp}.celltype{i_celltype}.scores_NSEM;
            
            scores_NSEM(find(scores_NSEM<-1)) = -1;
            
            if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
            if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
            plotstring_center = sprintf('%s.',colorstring);

            subplot(2,1,1); hold on;  
            plot(scores_WN(:,1), scores_WN(:,2), plotstring_center,  'markersize', MS_B);
            plot(scores_WN(:,1), scores_WN(:,2), plotstring_surround,  'markersize', MS_A);
            
            subplot(2,1,2); hold on
            plot(scores_NSEM(:,1), scores_NSEM(:,2), plotstring_center,  'markersize', MS_B);
            plot(scores_NSEM(:,1), scores_NSEM(:,2), plotstring_surround,  'markersize', MS_A);
        end

        subplot(2,1,1); hold on
        unity_line = linspace(0, 1,100);
        plot(unity_line,unity_line,'k')
        xlim([0,1]);
        ylim([0,1]);
        
        subplot(2,1,2); hold on
        unity_line = linspace(-1, 1,100);
        plot(unity_line,unity_line,'k')
        xlim([-1,1]);
        ylim([-1,1]);
    end

%orient landscape
%eval(sprintf('print -dpdf %s/ALLrast_vsGLM_BPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))




%{
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


%{
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
%}



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




%%
% RASTBPS_WRAP
% AKHEITMAN 2015-03-03: Considered Done
% WRAP of rastbps_comp
% Not STAND ALONE but calls "stand-alone" rastbps_comp
%
% GOAL:
% Compute the log-like of raster rate model with an added post spike term
% Normalization factor for the GLM which includes a conditioning term
% Normalizing with just a pure rate model not always sufficient for GLM
%
%
% INPUT:
% fittedGLM from NSEM_Home/GLM_Output_Analysis directory
% fittedGLM has the recorded raster,fitted ps filters, t_bin
%
% OUTPUT:
% base_crm_findPS_expA.mat file with structure "raster score"
% single structure with raster score for all cells in each experiment
% result put into NSEM_Home/Cell_Selection/Raster_Metrics directory

% AKHEITMAN 2015-03-09 
%    Added plot sections (compute, plotrasterscores, comparetoGLM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRECTORIES AND IDENTIFY WHICH CELLS TO COMPUTE

% RESET
clear; close all; clc

% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE
celltypes = [1 2]
exps = [1 2 3] 
debug = false;
compute = true; 
plotrasterscores = false;  
comparetoGLM = false;

comparetoGLM_withCONV = true;

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
 crm_type = 'base_crm_findPS_unitlikebasis';
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP LOAD COMPUTE SAVE
if compute
%%% REORGANIZE INTO A MORE APPROACHABLE STRUCTURE
if strcmp(crm_type, 'base_crm_findPS') || strcmp(crm_type, 'base_crm_findPS_HO')
    % HARD PS PARAMETERS  MATCHES GLMPARAMS
    t_bin0 = .00083275;
    ps_params.ms           = 100 ;     
    ps_params.filternumber = 20;
    ps_params.spacing      = pi/2;
    ps_params.bstretch     = .05;
    ps_params.alpha        = 0;
    ps_params.fratio       = .5;  
    ps_basis               = prep_spikefilterbasisGP(ps_params,t_bin0);
end


if strcmp(crm_type, 'base_crm_findPS_unitlikebasis')
    vec_1 = [1;1; zeros(38,1)];
    vec_2 = [1;1;1;1; zeros(76,1)];
    
    set_1 = NaN(40,20); pad_1 = zeros(80,20);
    set_2 = NaN(80,20); pad_2 = zeros(40,20);
    
    for i_rotate = 1:20
        rotate_1 = circshift(vec_1,2*(i_rotate-1) );
        rotate_2 = circshift(vec_2,4*(i_rotate-1) );
        
        set_1(:,i_rotate) =  rotate_1;
        set_2(:,i_rotate) =  rotate_2;
    end
    
    basis_1 = [set_1; pad_1];
    basis_2 = [pad_2; set_2];
    
    ps_basis = [basis_1 basis_2];
end


for i_exp = exps
    raster_scores = allcells{i_exp};
    raster_scores.stim_types = {'WN','NSEM'};
    raster_scores.celltypes  = {'ONP','OFFP'};
    
    raster_scores.notes.n1 = 'Both use 10 msec smoothing of raster spike trains';
    raster_scores.notes.n2 = 'UOP: unconditioned optimal performance';
    raster_scores.notes.n3 = 'CRM: conditioned rate model';
    raster_scores.notes.n4 = 'BPS: Bits Per Spike';
    raster_scores.notes.n5 = 'LPPS: logarithmic probability per second';
    raster_scores.notes.n6 = 'GLM: Generalized Linear Model';
    raster_scores.crm_type = crm_type;
    raster_scores.timestamp           = datestr(clock);
    raster_scores.mfile_name          = mfilename('fullpath');
    for i_celltype = celltypes
        % LOAD STIMULUS PARAMETERS / DEFINE CELL NUMBERS
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        map_type= 'mapPRJ';
        if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
        if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
      
        
        if exist('debug','var') && debug
            cellgroup = cellgroup(1:2);
        end
        
        scores_WN.uop_bps    = NaN(length(cellgroup),1);
        scores_WN.crm_bps    = NaN(length(cellgroup),1);
        scores_NSEM.uop_bps  = NaN(length(cellgroup),1);
        scores_NSEM.crm_bps  = NaN(length(cellgroup),1);
        scores_WN.uop_lpps   = NaN(length(cellgroup),1);
        scores_WN.crm_lpps   = NaN(length(cellgroup),1);
        scores_NSEM.uop_lpps = NaN(length(cellgroup),1);
        scores_NSEM.crm_lpps = NaN(length(cellgroup),1);
        
        scores_WN.glm_bps     = NaN(length(cellgroup),1);        
        scores_WN.glm_lpps    = NaN(length(cellgroup),1);
        scores_NSEM.glm_bps   = NaN(length(cellgroup),1);
        scores_NSEM.glm_lpps  = NaN(length(cellgroup),1);
        
        if strcmp(crm_type,'base_crm_findPS')
            scores_WN.psnote1   = 'Raster_Fitted_PSFilter';
            scores_WN.psnote2   = 'Exp of Filter is spike induced gain';
            scores_WN.ps_filter = cell(length(cellgroup),1);
            
            scores_NSEM.psnote1   = 'Raster_Fitted_PSFilter';
            scores_NSEM.psnote2   = 'Exp of Filter is spike induced gain';
            scores_NSEM.ps_filter = cell(length(cellgroup),1);
        end
        cellgroup
        for i_cell = 1:length(cellgroup)
            for i_stimtype = 1:2
                % LOAD FITTED GLM
                if i_stimtype == 1, stimtype = 'WN';   end
                if i_stimtype == 2, stimtype = 'NSEM'; end
                secondDir.exp_nm    = exp_nm;
                secondDir.stim_type = stimtype;
                secondDir.map_type  = 'mapPRJ';
                secondDir.fitname   = GLMType.fitname;
                Dirs.glmfitdir   = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
                clear fittedGLM cell_savename
                cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
                eval(sprintf('load %s/%s.mat', Dirs.glmfitdir, cell_savename))
                
                % PULL OUT COMPONENTS NECESSARY
                inp_rast = fittedGLM.xvalperformance.rasters.recorded;
                t_bin    = fittedGLM.t_bin;
                psfilter = fittedGLM.linearfilters.PostSpike.Filter; 
                
                glm_bps  = fittedGLM.xvalperformance.logprob_glm_bpspike;
                glm_lpps = fittedGLM.xvalperformance.logprob_glm_raw/(t_bin * size(inp_rast,2) ); 
                
                % COMPUTE BITS PER SPIKE FOR CONDITIONED RATE MODEL
                if strcmp(crm_type,'base_crm_importPS')
                    [uop_bits_perspike uop_logprob_persec crm_bits_perspike crm_logprob_persec pstar] = rastbps_comp_importPS(inp_rast,t_bin,psfilter);
                   
                elseif strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    [uop_bits_perspike uop_logprob_persec crm_bits_perspike crm_logprob_persec opt_params] = rastbps_comp_findPS(inp_rast,t_bin,ps_basis);
                     
                    display(sprintf('%s: Offset=%1.2e, ScaleRate=%1.2e', stimtype,opt_params.offset,opt_params.rate_drive));
                    
                elseif strcmp(crm_type,'base_crm_findPS_HO')
                    [uop_bits_perspike crm_bits_perspike pstar] = rastbps_comp_findPS(inp_rast,t_bin,options);
                   
                    %display(sprintf('%s Offset %1.2e,  Scale Rate %1.2e', stimtype,pstar(1),pstar(2)));
                end
                
                if i_stimtype == 1 
                    %WN_uop_bps = uop_bits_perspike;
                    %WN_crm_bps = crm_bits_perspike;
                    scores_WN.uop_bps(i_cell)  = uop_bits_perspike;
                    scores_WN.crm_bps(i_cell)  = crm_bits_perspike;
                    scores_WN.uop_lpps(i_cell) = uop_logprob_persec;
                    scores_WN.crm_lpps(i_cell) = crm_logprob_persec;   
                    scores_WN.glm_bps(i_cell)  = glm_bps;
                    scores_WN.glm_lpps(i_cell) = glm_lpps;   
                elseif i_stimtype == 2
                    %NSEM_uop_bps = uop_bits_perspike;
                    %NSEM_crm_bps = crm_bits_perspike;
                    scores_NSEM.uop_bps(i_cell)  = uop_bits_perspike;
                    scores_NSEM.crm_bps(i_cell)  = crm_bits_perspike;
                    scores_NSEM.uop_lpps(i_cell) = uop_logprob_persec;
                    scores_NSEM.crm_lpps(i_cell) = crm_logprob_persec;  
                    scores_NSEM.glm_bps(i_cell)  = glm_bps;
                    scores_NSEM.glm_lpps(i_cell) = glm_lpps;   
                end
                
                if strcmp(crm_type,'base_crm_findPS') || strcmp(crm_type ,'base_crm_findPS_unitlikebasis')
                    if i_stimtype == 1
                        scores_WN.ps_filter{i_cell}   = opt_params.ps_filter;
                    elseif i_stimtype == 2
                        scores_NSEM.ps_filter{i_cell} = opt_params.ps_filter;
                    end
                end
                
                
            end
            display(sprintf('### %s %s: CRM: LogProb of WN-NSEM: %1.2e ###',...
                expname,cell_savename, (scores_WN.crm_lpps(i_cell)-scores_NSEM.crm_lpps(i_cell))  ));
            %display(sprintf('### %s %s: UOP: LogProb of WN-NSEM: %1.2e ###',...
            %    expname,cell_savename,  (scores_WN.uop_lpps(i_cell)-scores_NSEM.uop_lpps(i_cell)) ));
            
        end
        raster_scores.celltype{i_celltype}.scores_WN = scores_WN;
        raster_scores.celltype{i_celltype}.scores_NSEM = scores_NSEM;
       
        if exist('debug','var') && debug
            eval(sprintf('save %s/%s_%s_DBUG.mat raster_scores',savedir,crm_type,exp_nm));
        else
            eval(sprintf('save %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
        end
    end
end

end

%% PLOTTING RAW RASTER SCORES
if plotrasterscores
% Plot WN-NSEM comparison for these base predictability metrics
colors = {'r','g','b','c'};
MS_A = 12;
MS_B = 16;

for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
    colorstring = colors{i_exp};
    
    for i_celltype = celltypes
        figure; clf;
        if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end
        
        scores_WN   = raster_scores.celltype{i_celltype}.scores_WN;
        scores_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM;
        
        if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
        if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
        plotstring_center = sprintf('%s.',colorstring);
        
        subplot(5,1,1)
        set(gca, 'fontsize', 8); axis off
        c = 0; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s, Conditioning Type: %s',exp_nm, celltype,crm_type),'interpreter','none');
        c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
        if strcmp(crm_type, 'base_crm_importPS')
            c=c+1; text(-.1, 1-0.15*c,sprintf('PS filter imported from GLM Fit: %s', GLMType.fitname), 'interpreter','none');
        end
        c=c+1; text(-.1, 1-0.15*c,sprintf('Rasters Ability to Predict Itself: Pure rate model and Conditioned Rate Model'));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Logarithmic Prob Per Sec - Raw attempt at probability of raster given model,more negative is worse'));
        c=c+1; text(-.1, 1-0.15*c,sprintf('BPS: Bits Per Spike - Improvement from Flat Rate Model, more positive is better'));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
    
        for i_metric = 1:4
            
            if i_metric == 1
                vals_WN   = scores_WN.uop_lpps; 
                vals_NSEM = scores_NSEM.uop_lpps;
                met_name  = 'Unconditioned_LogProbperSec';
                plot_index = [3,5];
            end
            if i_metric == 2
                vals_WN   = scores_WN.crm_lpps; 
                vals_NSEM = scores_NSEM.crm_lpps;
                met_name  = 'Conditioned_LogProbperSec';
                plot_index = [4,6];
            end
            if i_metric == 3
                vals_WN   = scores_WN.uop_bps; 
                vals_NSEM = scores_NSEM.uop_bps;
                met_name  = 'Unconditioned_BitsPerSpike';
                plot_index = [7,9];
            end
            
            if i_metric == 4
                vals_WN   = scores_WN.crm_bps; 
                vals_NSEM = scores_NSEM.crm_bps;
                met_name  = 'Conditioned_BitsPerSpike';
                plot_index = [8,10];
            end
            subplot(5,2,plot_index); hold on; set(gca,'fontsize',8)
            maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',met_name),'interpreter','none')

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
           % xlabel('White Noise');
           % ylabel('Natural Scenes');
        end
        
        orient landscape
        eval(sprintf('print -dpdf %s/rast_selfprediction_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))
    end
end


end
%% COMPARE RASTER SCORES TO GLM PREDICTABILITY

if comparetoGLM
colors = {'r','g','b','c'};
MS_A = 12;
MS_B = 16;

for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
    colorstring = colors{i_exp};
    
    for i_celltype = celltypes
        for i_mettype = 1:2
            clf;
            if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
            if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end

            scores_WN   = raster_scores.celltype{i_celltype}.scores_WN;
            scores_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM;

            if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
            if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
            plotstring_center = sprintf('%s.',colorstring);

            subplot(5,1,1)
            set(gca, 'fontsize', 10); axis off
            c = 0; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s, Conditioning Type: %s',exp_nm, celltype,crm_type),'interpreter','none');
            c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
            c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
            c=c+1; text(-.1, 1-0.15*c,sprintf('GLM Type: %s', GLMType.fitname),'interpreter','none');
            if i_mettype == 1
            c=c+1; text(-.1, 1-0.15*c,sprintf('BPS: Bits Per Spike, Improvement from Flat Rate Model: Relative Performance, Postive is Better'));
            elseif i_mettype == 2
                c=c+1; text(-.1, 1-0.15*c,sprintf('LPPS: Logarthimic Probability Per Second. Pure performance, less negative is better'));
            end
            c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
            c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;

            for i_metric = 1:4
                if i_metric == 1, met_name = 'Unconditioned Subtract'; plot_index = [3,5]; end
                if i_metric == 2, met_name = 'Unconditioned Divide';   plot_index = [4,6]; end
                if i_metric == 3, met_name = 'Conditioned Subtract'; plot_index = [7,9]; end
                if i_metric == 4, met_name = 'Conditioned Divide';   plot_index = [8,10]; end

                % GLM performance
                if i_mettype == 1
                    glm_WN   = scores_WN.glm_bps; 
                    glm_NSEM = scores_NSEM.glm_bps;   
                elseif i_mettype == 2
                    glm_WN   = scores_WN.glm_lpps; 
                    glm_NSEM = scores_NSEM.glm_lpps;
                end

                % Define Normalizing Vectors
                if i_mettype == 1
                    if i_metric == 1 || i_metric == 2;
                        norm_WN   = scores_WN.uop_bps; 
                        norm_NSEM = scores_NSEM.uop_bps;                
                    end
                    if i_metric == 3 || i_metric == 4;
                        norm_WN   = scores_WN.crm_bps; 
                        norm_NSEM = scores_NSEM.crm_bps;                
                    end
                elseif i_mettype == 2
                    if i_metric == 1 || i_metric == 2;
                        norm_WN   = scores_WN.uop_lpps; 
                        norm_NSEM = scores_NSEM.uop_lpps;                
                    end
                    if i_metric == 3 || i_metric == 4;
                        norm_WN   = scores_WN.crm_lpps; 
                        norm_NSEM = scores_NSEM.crm_lpps;                
                    end
                end

                % Construct Scores
                if i_metric == 1 || i_metric == 3
                    vals_WN    = glm_WN - norm_WN;
                    vals_NSEM  = glm_NSEM - norm_NSEM;
                elseif i_metric == 2 || i_metric == 4
                    vals_WN    = (glm_WN)./norm_WN;
                    vals_NSEM  = (glm_NSEM)./norm_NSEM;
                end

                % Throw out Garbage
                bad_ind = union(find(abs(vals_WN) == Inf),find(abs(vals_NSEM) == Inf));
                vals_WN(bad_ind) = [];
                vals_NSEM(bad_ind) = [];

                %  
                subplot(5,2,plot_index); hold on; set(gca,'fontsize',10)
                maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
                minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
                plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
                plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
                title(sprintf('%s',met_name),'interpreter','none')

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
                %xlabel('White Noise');
                %ylabel('Natural Scenes');
            end

            orient landscape
            if i_mettype == 1
                eval(sprintf('print -dpdf %s/rast_vsGLM_BPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))
            elseif i_mettype == 2
                eval(sprintf('print -dpdf %s/rast_vsGLM_LPPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))
            end
        end
    end
end

end
%%


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



%% OLD WAY
%{
if comparetoGLM
colors = {'r','g','b','c'};
MS_A = 12;
MS_B = 16;

for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,crm_type,exp_nm)); 
    colorstring = colors{i_exp};
    
    for i_celltype = celltypes
        clf;
        if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end
        
        scores_WN   = raster_scores.celltype{i_celltype}.scores_WN;
        scores_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM;
       
        if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
        if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
        plotstring_center = sprintf('%s.',colorstring);
        
        subplot(5,1,1)
        set(gca, 'fontsize', 10); axis off
        c = 0;
        c=c+1; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s, Conditioning Type: %s',exp_nm, celltype,crm_type),'interpreter','none');
        c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
        c=c+1; text(-.1, 1-0.1*c,sprintf('GLM Performance Normalized by Raster Scores: Pure rate model and Conditioned Rate Model'));
        c=c+1; text(-.1, 1-0.1*c,sprintf('BPS: Bits Per Spike - Improvement from Flat Rate Model, more positive is better'));
        c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
        c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
    
        for i_metric = 1:4
            if i_metric == 1, met_name = 'Unconditioned Subtract'; plot_index = [3,5]; end
            if i_metric == 2, met_name = 'Unconditioned Divide';   plot_index = [4,6]; end
            if i_metric == 3, met_name = 'Conditioned Subtract'; plot_index = [7,9]; end
            if i_metric == 4, met_name = 'Conditioned Divide';   plot_index = [8,10]; end
            
            % GLM performance
            glm_WN   = scores_WN.glm_bps; 
            glm_NSEM = scores_NSEM.glm_bps;   
            
            % Define Normalizing Vectors
            if i_metric == 1 || i_metric == 2;
                norm_WN   = scores_WN.uop_bps; 
                norm_NSEM = scores_NSEM.uop_bps;                
            end
            if i_metric == 3 || i_metric == 4;
                norm_WN   = scores_WN.crm_bps; 
                norm_NSEM = scores_NSEM.crm_bps;                
            end
            
            % Construct Scores
            if i_metric == 1 || i_metric == 3
                vals_WN    = glm_WN - norm_WN;
                vals_NSEM  = glm_NSEM - norm_NSEM;
            elseif i_metric == 2 || i_metric == 4
                vals_WN    = (glm_WN)./norm_WN;
                vals_NSEM  = (glm_NSEM)./norm_NSEM;
            end
            
            % Throw out Garbage
            bad_ind = union(find(abs(vals_WN) == Inf),find(abs(vals_NSEM) == Inf));
            vals_WN(bad_ind) = [];
            vals_NSEM(bad_ind) = [];
            
            %  
            subplot(5,2,plot_index); hold on; set(gca,'fontsize',10)
            maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',met_name),'interpreter','none')

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
        eval(sprintf('print -dpdf %s/rast_vsGLM_BPS_%s_%s_%s.pdf', savedir, crm_type, exp_nm, celltype))
    end
end
%}

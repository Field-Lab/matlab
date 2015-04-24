%%% PURPOSE %%%
% Collect and plot GLM_Scores in an easy and managebale fashion

%%% COMPUTATIONS %%%
% NONE

% AKHEITMAN 2014-12-09 start

%metric_type = 'Viktor_Spike';
%metric_type = 'bits_per_spike';
clear; clc; close all
exps = 1:4; stimtypes = 1:2; celltypes = 1:2;
i_exp = 1; i_stimtype = 1; i_celltype = 1; i_cell =1;


%metric_type = 'bits_per_spike';
%metric_type = 'fracvar_avgsignals';
metric_type = 'Viktor_Spike'
if strcmp(metric_type,'fracvar_avgsignals')
    normalization = false;
end
if strcmp(metric_type,'Viktor_Spike')
    normalization = true;
end

%if strcmp(metric_type,'bits_per_spike_opttime')
%    normalization = true;
%end


BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
glm_scores = allcells;


% HARD PARAMETERS OF THE RASTER WE ALREADY HAVE SET
hard_params.raster_params.bindur         = .00083275;
hard_params.raster_params.bins_per_frame = 10;
hard_params.map_type = 'mapPRJ';

% SPECIFY GLMTYPE TO GET A FITNAME
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.nullpoint = 'mean';  GLMType.debug = false;
GLMType.CONVEX = true; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.specialchange = false;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.fitname  = GLM_fitname(GLMType); 
fitname = GLMType.fitname;

% CREATE SAVEDIR
Dirs.glmscores = sprintf('%s/GLM_Performance/%s', BD.Cell_Selection,fitname);
Dirs.savedir   = sprintf('%s/%s', Dirs.glmscores, metric_type);
if ~exist(Dirs.savedir, 'dir'), mkdir(Dirs.savedir); end


filename = sprintf('%s/%s',Dirs.savedir, metric_type)


%% SAVE STRUCTURE IF NECESSARY
if ~exist(sprintf('%s.mat',filename), 'file')
for i_exp = exps
    for i_stimtype = stimtypes   
        for i_celltype = celltypes
            %% Experiment Dependent Parameters
            
            % CLEAN UP
            clear StimPars raster_params secondDir
            % LOAD STIMULUS PARAMETERS / DEFINE CELL NUMBERS
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            map_type= hard_params.map_type;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            
            % FITTED GLM
            secondDir.exp_nm = exp_nm;
            secondDir.stim_type = stimtype;
            secondDir.map_type  = 'mapPRJ';
            secondDir.fitname   = fitname;
            Dirs.fittedGLM      = NSEM_secondaryDirectories('savedir_GLMfit', secondDir);  clear secondDir
            Dirs.crossval_load  = sprintf('%s/crossval_%s', Dirs.fittedGLM, metric_type);
            
            
            % INITIALIZE THE SOLUTION STRUCTURE
            cid = cellgroup(1);
            cell_savename = sprintf('%s_%d', celltype,cid);
            eval(sprintf('load %s/crossvalperf_%s.mat', Dirs.crossval_load, cell_savename));
            scores.metric_type = metric_type;
            scores.metric_note = 'rows are cells, columns are different time scales';
            scores.cells       = cellgroup;
            if strcmp(metric_type,'fracvar_avgsignals')
                scores.binscales   = crossval_perf.scores.sigma_bin';
            elseif strcmp(metric_type, 'bits_per_spike')
                scores.binscales   = crossval_perf.scores.sigma_bin_norm;
            elseif strcmp(metric_type, 'Viktor_Spike')
                scores.binscales   = crossval_perf.scores.Viktor_Time_Bins;
            end
            scores.values_raw = zeros(length(scores.binscales), length(cellgroup) );
            if normalization
                scores.values_normsubtract = zeros(length(scores.binscales), length(cellgroup) );
                scores.values_normdivide   = zeros(length(scores.binscales), length(cellgroup) );
            end
            
            %%
            for i_cell = 1:length(cellgroup)
                % LOAD STORED RASTER METRICS
                cid = cellgroup(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
                display(sprintf('Working on Exp: %s Stim: %s Cell: %s', exp_nm,stimtype, cell_savename));               
                eval(sprintf('load %s/crossvalperf_%s.mat', Dirs.crossval_load, cell_savename));
                
                scores.values_raw(:,i_cell) =  crossval_perf.scores.metric_raw;
                if normalization
                    scores.values_normsubtract(:,i_cell) = crossval_perf.scores.metric_normsubtract;
                    scores.values_normdivide(:,i_cell)   = crossval_perf.scores.metric_normdivide;
                end                
            end
            glm_scores{i_exp}.stim_type{i_stimtype}.celltype{i_celltype}.scores = scores; 
        end
    end
end
eval(sprintf('save %s.mat glm_scores',filename));
end

%%




eval(sprintf('load %s.mat glm_scores',filename));
eval(sprintf('load %s/GoodSim_Rules.mat', BD.Cell_Selection))
binscales = glm_scores{1}.stim_type{1}.celltype{1}.scores.binscales;
norm_type = 'normsubtract';
norm_type = 'raw';
%norm_type = 'normdivide';
%rule = 'allcells';
rule = 'glmconverge_rule3';

for i_bin = 1:length(binscales)
    clf;
    subplot(5,1,1);
    MS = 10;
    axis off
    c = 0;
    text(-.1, 1,sprintf('Xaxis: WN Values,  Yaxis: NSEM Values,   ALL CELLS' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Fit Type: %s', fitname),'interpreter','none');
    c=c+1; text(-.1, 1-0.1*c,'Color are experiments, dots ONP, asterisk OFFP');
    c=c+1; text(-.1, 1-0.1*c,'0 value means worse than steady firing rate,1 means unconditioned optimum');
    for i_exp = 1:4
        if i_exp == 1; basecolor = 'r'; end
        if i_exp == 2; basecolor = 'g'; end
        if i_exp == 3; basecolor = 'b'; end
        if i_exp == 4; basecolor = 'c'; end
        
        for i_celltype = 1:2
            if i_celltype == 1, marktype  = '.'; cellgroup = glm_scores{i_exp}.ONP; end
            if i_celltype == 2, marktype  = '*'; cellgroup = glm_scores{i_exp}.OFFP;end
            
            if strcmp(norm_type, 'raw')
                WNvals      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw(i_bin,:);
                NSEMvals    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw(i_bin,:);
            elseif strcmp(norm_type, 'normdivide')
                WNvals      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide(i_bin,:);
                NSEMvals    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide(i_bin,:);
            elseif strcmp(norm_type, 'normsubtract')
                WNvals      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract(i_bin,:);
                NSEMvals    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract(i_bin,:);
            end
            
            
            if strcmp(rule, 'glmconverge_rule3')
                [dummy out_index] = intersect(cellgroup , union(Rules{3}.byexpnm{i_exp}.outcells_NSEM , Rules{3}.byexpnm{i_exp}.outcells_WN) )
            end
                
            
            badind  = union(find(NSEMvals ==  NaN),union(find(NSEMvals == -Inf),find(NSEMvals ==  Inf)));
            
            %{
            WNvals(badind) = [];
            NSEMvals(badind) = [];
            %}
            
            
            if strcmp(metric_type, 'bits_per_spike')
                WNvals  ( WNvals   <=  -1) = -1;
                NSEMvals( NSEMvals <=  -1) = -1;
            end
         
            
            subplot(5,2, (i_exp*2 + i_celltype))
            hold on;
            max_val = max(max(WNvals),max(NSEMvals));
            min_val = min(min(WNvals),min(NSEMvals));
            

            %set(gca,'xlim',[min_val max_val]); set(gca,'ylim',[min_val max_val]);
            %set(gca,'xlim',[-1 1]); set(gca,'ylim',[-1 1]);
            
            %set(gca,'xlim',[-1 max_val]); set(gca,'ylim',[-1 max_val]);
            set(gca,'xlim',[min_val max_val]); set(gca,'ylim',[min_val max_val]);
            %set(gca,'xlim',[-1 0]); set(gca,'ylim',[-1 1]);
            plotstring = sprintf('%s%s',basecolor,marktype);
            plot(WNvals, NSEMvals, plotstring,'markersize',MS);
            plot(WNvals(out_index), NSEMvals(out_index),'k.','markersize',MS)
            plot(linspace(min_val,max_val,100), linspace(min_val,max_val,100), 'k')
        end
    end

    orient tall
    eval(sprintf('print -dpdf %s/WN_V_NSEM_%s_%s_Bin%d.pdf', Dirs.savedir,norm_type, rule,binscales(i_bin)))
    
    
    
end
%%
% Pick Maximum Value FOr Subtraction %
for i_bin = 1
    clf;
    subplot(5,1,1);
    MS = 10;
    axis off
    c = 0;
    text(-.1, 1,sprintf('Xaxis: WN Values,  Yaxis: NSEM Values,   ALL CELLS' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Fit Type: %s', fitname),'interpreter','none');
    c=c+1; text(-.1, 1-0.1*c,'Color are experiments, dots ONP, asterisk OFFP');
    c=c+1; text(-.1, 1-0.1*c,'0 value means worse than steady firing rate,1 means unconditioned optimum');
    for i_exp = 1:4
        if i_exp == 1; basecolor = 'r'; end
        if i_exp == 2; basecolor = 'g'; end
        if i_exp == 3; basecolor = 'b'; end
        if i_exp == 4; basecolor = 'c'; end
        
        for i_celltype = 1:2
            if i_celltype == 1, marktype  = '.'; end
            if i_celltype == 2, marktype  = '*'; end
            
            
            WNvals_normfactor   = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw - glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract;
            NSEMvals_normfactor = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw - glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract;
            if strcmp(norm_type, 'raw')
                WNvals_full      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_raw;
                NSEMvals_full    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_raw;
            elseif strcmp(norm_type, 'normdivide')
                WNvals_full      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normdivide;
                NSEMvals_full    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normdivide;
            elseif strcmp(norm_type, 'normsubtract')
                WNvals_full      = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values_normsubtract;
                NSEMvals_full    = glm_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values_normsubtract;
            end
            
            cellgroup = glm_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.cells;
            WNvals   = zeros(1,length(cellgroup));
            NSEMvals = zeros(1,length(cellgroup));
            
            for i_cell = 1:length(cellgroup);
                [dummy,WN_index]   = max(WNvals_normfactor(:,i_cell));
                [dummy,NSEM_index] = max(NSEMvals_normfactor(:,i_cell));
                WNvals(i_cell)   = WNvals_full(WN_index,i_cell);
                NSEMvals(i_cell) = NSEMvals_full(NSEM_index,i_cell);
            end
            
            badind         = union(find(NSEMvals ==  NaN),union(find(NSEMvals == -Inf),find(NSEMvals ==  Inf)));
            WNvals(badind) = [];
            NSEMvals(badind) = [];
            
            subplot(5,2, (i_exp*2 + i_celltype))
            hold on;
            
            WNvals  ( WNvals <=-1)   = -1;
            NSEMvals( NSEMvals <=-1) = -1;
            
            max_val = max(max(WNvals),max(NSEMvals));
            min_val = min(min(WNvals),min(NSEMvals));
            

            %set(gca,'xlim',[min_val max_val]); set(gca,'ylim',[min_val max_val]);
            %set(gca,'xlim',[-1 1]); set(gca,'ylim',[-1 1]);
            
            %set(gca,'xlim',[-1 max_val]); set(gca,'ylim',[-1 max_val]);
            set(gca,'xlim',[min_val max_val]); set(gca,'ylim',[min_val max_val]);
            %set(gca,'xlim',[-1 0]); set(gca,'ylim',[-1 1]);
            plotstring = sprintf('%s%s',basecolor,marktype);
            plot(WNvals, NSEMvals, plotstring,'markersize',MS);
            plot(linspace(min_val,max_val,100), linspace(min_val,max_val,100), 'k')
        end
    end

    orient tall
    eval(sprintf('print -dpdf %s/WN_V_NSEM_%s_%s_MAXTIME.pdf', Dirs.savedir,norm_type,rule))
    
    
    
end


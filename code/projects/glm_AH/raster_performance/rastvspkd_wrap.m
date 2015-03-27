%%
% VSPKD_WRAP
% AKHEITMAN 2015-03-09

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIRECTORIES AND IDENTIFY WHICH CELLS TO COMPUTE

% RESET
clear; close all; clc

% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE
celltypes    = [1 2]
exps         = [1 2 3 4] %3 4] 
debug        = false;
compute      = false; 
plotrasterscores = true;  
%comparetoGLM = true;

% SET DIRECTORIES / USING GLM becaause it has the raster.
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
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

% Victor Types
pars.msec    =  50; 


met_type = sprintf('vspkd_msec%d', pars.msec);
if exist('debug','var') && debug
    met_type = sprintf('%s_DBUG',met_type);
end

% Victor Bin 
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP LOAD COMPUTE SAVE
if compute
    
%%% REORGANIZE INTO A MORE APPROACHABLE STRUCTURE
for i_exp = exps
    raster_scores = allcells{i_exp};
    raster_scores.metric_type= met_type;
    raster_scores.stim_types = {'WN','NSEM'};
    raster_scores.celltypes  = {'ONP','OFFP'};
    raster_scores.notes.n1   = 'Using Victor Spike Code from http://www-users.med.cornell.edu/~jdvicto/spkdm.html';
    raster_scores.timestamp  = datestr(clock);
    raster_scores.mfile_name = mfilename('fullpath');
    %%
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
        scores_WN.rast_vspkd      = NaN(length(cellgroup),1);
        scores_NSEM.rast_vspkd    = NaN(length(cellgroup),1);
        %%
        for i_cell = 1:length(cellgroup)
            tStart = tic;
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
                rec_rast = fittedGLM.xvalperformance.rasters.recorded;
                t_bin    = fittedGLM.t_bin;
                totalreps  = size(rec_rast,1);
                reps       = floor(totalreps/2);
                subsets    = randsample((1:totalreps), 2*reps);
                a = intersect(subsets(1:reps), subsets((reps+1):2*reps));
                if ~isempty(a) 
                    error('subrasters overlap'); 
                end
                clear a
                subrast_1   = rec_rast(subsets(1:reps),:);
                subrast_2   = rec_rast(subsets((reps+1):2*reps),:);
                
                if exist('debug','var') && debug
                    subrast_1 = subrast_1(1:3,:);
                    subrast_2 = subrast_2(1:3,:);
                end
                
                
                timescale   = pars.msec / 1000;
                [vspkd_mean vspkd_std] = rastvspkd_comp(subrast_1,subrast_2,timescale,t_bin);

                % COMPUTE VICTOR SPIKE 
                if i_stimtype == 1 
                    scores_WN.rast_vspkd(i_cell)      = vspkd_mean;    
                elseif i_stimtype == 2
                    scores_NSEM.rast_vspkd(i_cell)    = vspkd_mean;  
                end                
            end
            duration = toc(tStart);
            display(sprintf('%s %s, WN: %1.2e, NSEM: %1.2e, Comp Time in Minutes: %1.2e',...
                expname,cell_savename, scores_WN.rast_vspkd(i_cell),scores_NSEM.rast_vspkd(i_cell), duration/60  ));
         
        end
        raster_scores.celltype{i_celltype}.scores_WN = scores_WN;
        raster_scores.celltype{i_celltype}.scores_NSEM = scores_NSEM;
       
       
        eval(sprintf('save %s/%s_%s.mat raster_scores',savedir,met_type,exp_nm)); 

    end
end

end


%%

if plotrasterscores
% Plot WN-NSEM comparison for these base predictability metrics
colors = {'r','g','b','c'};
MS_A = 12;
MS_B = 16;

for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,met_type,exp_nm)); 
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
        set(gca, 'fontsize', 8); axis off
        c = 0; text(-.1, 1-0.15*c,sprintf('Experiment %s:  Celltype: %s',exp_nm, celltype),'interpreter','none');
        c=c+1; text(-.1, 1-0.15*c,sprintf('Victor Spike divided by Spike Count: %s', met_type),'interpreter','none');
        c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
        c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
    

            subplot(5,1,[2 3 4 5]); hold on; set(gca,'fontsize',8)
            vals_WN = scores_WN.rast_vspkd;
            vals_NSEM = scores_NSEM.rast_vspkd;
            
            maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',met_type),'interpreter','none')

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

        
        orient landscape
        eval(sprintf('print -dpdf %s/%s_PLOT_%s_%s.pdf', savedir, met_type, exp_nm, celltype))
    end
end


clf;
subplot(5,1,1)
set(gca, 'fontsize', 8); axis off
c = 0; text(-.1, 1-0.15*c,sprintf('Victor Spike divided by Spike Count: %s', met_type),'interpreter','none');
c=c+1; text(-.1, 1-0.15*c,sprintf('X-Axis: White Noise Values, Y-Axis: NSEM Values'));
c=c+1; text(-.1, 1-0.15*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-0.15*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none') ;
subplot(5,1,[2 3 4 5]); hold on; set(gca,'fontsize',8)
for i_exp = exps
    clear raster_scores
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    eval(sprintf('load  %s/%s_%s.mat raster_scores',savedir,met_type,exp_nm)); 
    colorstring = colors{i_exp};
    
    for i_celltype = celltypes
        if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end
        
        scores_WN   = raster_scores.celltype{i_celltype}.scores_WN;
        scores_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM;
        
        if i_celltype == 1, plotstring_surround = sprintf('%so', colorstring);  end
        if i_celltype == 2, plotstring_surround = sprintf('%ss', colorstring);  end
        plotstring_center = sprintf('%s.',colorstring);
        
        
    

            
            vals_WN = scores_WN.rast_vspkd;
            vals_NSEM = scores_NSEM.rast_vspkd;
            
            maxval = max(max(vals_WN(:)), max(vals_NSEM(:)));
            minval = min(min(vals_WN(:)), min(vals_NSEM(:)));
            plot(vals_WN, vals_NSEM, plotstring_center,  'markersize', MS_B);
            plot(vals_WN, vals_NSEM, plotstring_surround,'markersize', MS_A);
            title(sprintf('%s',met_type),'interpreter','none')

            if exist('label_points','var') && label_points
                for i_cell = 1:length(cids)
                    lab = sprintf('%d', cids(i_cell));
                    text(vals(1,i_cell),vals(2,i_cell),lab);
                end
            end
            

        
        
    end
end
unity_line = linspace(0,2,100);
plot(unity_line,unity_line,'k')
xlim([0,2]);
ylim([0,2]);
xlabel('White Noise');
ylabel('Natural Scenes');
orient landscape
eval(sprintf('print -dpdf %s/%s_PLOTALL.pdf', savedir, met_type))

end

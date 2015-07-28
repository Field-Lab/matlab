% AKHEITMAN 2015-07-16
% Define convergence without ignoring improvement in first half
%{
clear; close all; clc
exps = [1 2 3 4];
glm_settings{1}.type = 'filter_mode';
glm_settings{1}.name = 'rk1';
special_arg = {};
runoptions.compute_50fit     = false;
runoptions.compute_halfratio = true;
runoptions.plot = true;
shortname = 'halfratio_standardrk1';
eval_conv_fit50(exps,glm_settings,special_arg,shortname,runoptions);
%}


function eval_conv_fit50(exps,glm_settings,special_arg,shortname,runoptions);
%% Bookkeep
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
GLMType = GLM_settings('default',glm_settings);
GLMType.fitname    = GLM_fitname(GLMType);
if exist('special_arg','var')
    for i_arg = 1:length(special_arg);
        GLMType.fitname = sprintf('%s/%s', GLMType.fitname, special_arg{i_arg});
    end
end

fit50_GLMType   = GLMType;
fit50_GLMType.fitname = sprintf('%s/%s', fit50_GLMType.fitname, 'fit_50percent');


fit0_glm_settings{1}.type = 'filter_mode';
fit0_glm_settings{1}.name = 'nostim';
fit0_glm_settings{2}.type = 'PostSpikeFilter';
fit0_glm_settings{2}.name =  'OFF';
fit0_GLMType = GLM_settings('default',fit0_glm_settings);
fit0_GLMType.fitname    = GLM_fitname(fit0_GLMType);

Dirs.convergence = sprintf('%s/%s/convergence', BD.GLM_output_analysis,GLMType.fitname);
if ~exist(Dirs.convergence, 'dir'), mkdir(Dirs.convergence); end

%% Compute 50 percent score
if runoptions.compute_50fit    
metric_type.name  = 'objval';
metric_type.note = 'Computed Objective Value from parameter fitting routine'    
model.fitname = fit50_GLMType.fitname;
savedir = sprintf('%s/%s', BD.GLM_output_analysis, model.fitname)

for i_exp = exps
    exp_nm  = allcells{i_exp}.exp_nm;
    filename = sprintf('%s/%s_%s.mat', savedir,'objval', exp_nm)
    expname = allcells{i_exp}.expname;
    
    
    % Create the aggregated scores structure
    if ~exist(filename)
        aggregated_scores = allcells{i_exp};  
        aggregated_scores.metric = metric_type.name;
        aggregated_scores.metric_note = metric_type.note;
        for i_celltype = [1,2]
            if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
            if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
            aggregated_scores.celltype{i_celltype}.scores_WN   = NaN(length(cellgroup),1);
            aggregated_scores.celltype{i_celltype}.scores_NSEM = NaN(length(cellgroup),1);
        end   
    else
        display(sprintf('checking for updates to %s on %s', model.fitname,exp_nm));
        eval(sprintf('load %s/%s_%s.mat aggregated_scores', savedir,metric_type.name, exp_nm));
        combined_scores = [aggregated_scores.celltype{1}.scores_WN; aggregated_scores.celltype{1}.scores_NSEM;...
           aggregated_scores.celltype{2}.scores_WN; aggregated_scores.celltype{1}.scores_NSEM];
        exist_nan = sum(isnan(combined_scores));
        if exist_nan == 0
           display('No Updates Necessary .. All scores Present')
           continue
        end
    end
    
    for i_celltype = [1,2]
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        
        for i_stimtype = 1:2
            % Pull the scores
            if i_stimtype == 1, rawscores = aggregated_scores.celltype{i_celltype}.scores_WN; stimtype = 'WN';  end
            if i_stimtype == 2, rawscores = aggregated_scores.celltype{i_celltype}.scores_NSEM; stimtype = 'NSEM'; end
            
            % Find the correct load directory
            secondDir.exp_nm    = exp_nm;
            secondDir.fitname = model.fitname;
            secondDir.map_type  = 'mapPRJ';
            secondDir.stim_type = stimtype;
            if i_stimtype == 1, secondDir.stim_type = 'WN'; end
            if i_stimtype == 2, secondDir.stim_type = 'NSEM'; end  
            loaddir = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);           
            organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
            clear secondDir

            GLMType.fit_type = stimtype;

            % Load and process stimulus
            [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
            [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'fitmovie');
            GLMType.fitmoviefile  = origmatfile;        
            fitmovie  = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv); 
            
            
            stillNAN = find(isnan(rawscores));    
            % Run through cells that still have NaN as the entry
            for i_cell = 1:length(stillNAN);
                i_index = stillNAN(i_cell);
                clear crossval_rawmetrics cell_savename
                cid = cellgroup(i_index);
                cell_savename = sprintf('%s_%d', celltype,cid);
                matfilename = sprintf('%s/%s.mat', loaddir, cell_savename);
                
                % Load Blocked-Spikes and make fitspikes_concat
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', organizedspikesdir, cell_savename));
                fitspikes.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                
                
                % Compute metric if necessary
                 if exist(matfilename)
                    display(sprintf('LOADING %s %s', expname,cell_savename));
                    eval(sprintf('load %s', matfilename));
                    rawscores(i_index);
                    t_bin = fittedGLM.t_bin;
                    bins  = fittedGLM.bins_per_frame * size(fitmovie,3);
                    pstar = fittedGLM.rawfit.opt_params;
                    [lcif_fit.nonstim] = subR_lcif_nonstim(pstar, fittedGLM.GLMType,fittedGLM.GLMPars,fitspikes,t_bin,bins);
                    [lcif_fit.stim.preNL objval] = subR_findobj_lcifstim(lcif_fit.nonstim.total, pstar, ...
                        fittedGLM.GLMType, fittedGLM.GLMPars, fitspikes, fitmovie, inputstats, fittedGLM.cellinfo, t_bin,bins);  
                    display(sprintf('objval is %d', objval))
                    if strcmp(metric_type.name,'objval')
                        rawscores(i_index) = objval;
                    end
                 end
            end

            if i_stimtype == 1
                aggregated_scores.celltype{i_celltype}.scores_WN = rawscores;
            end
            if i_stimtype == 2
                aggregated_scores.celltype{i_celltype}.scores_NSEM = rawscores;
            end
        end
    end
    aggregated_scores.timestamp       = datestr(clock);
    aggregated_scores.code_name       = mfilename('fullpath');
    eval(sprintf('save %s/%s_%s.mat aggregated_scores', savedir,metric_type.name, exp_nm))
end
end

%% Compute Half Ratio

if runoptions.compute_halfratio
agg_scores = allcells;

Dirs.load_fit0    = sprintf('%s/%s', BD.GLM_output_analysis,  fit0_GLMType.fitname);
Dirs.load_fit50   = sprintf('%s/%s', BD.GLM_output_analysis, fit50_GLMType.fitname);
Dirs.load_fitFULL = sprintf('%s/%s', BD.GLM_output_analysis,       GLMType.fitname);


for i_exp = exps
    % load data structures
    exp_nm = allcells{i_exp}.exp_nm;    
    eval(sprintf('load %s/objval_%s', Dirs.load_fit0, exp_nm));
    objval.fit0 = aggregated_scores; clear aggregated_scores 
    eval(sprintf('load %s/objval_%s', Dirs.load_fit50, exp_nm));
    objval.fit50 = aggregated_scores; clear aggregated_scores 
    eval(sprintf('load %s/objval_%s', Dirs.load_fitFULL, exp_nm));
    objval.fitFULL = aggregated_scores; clear aggregated_scores 
    
    % initialize structure
    agg_scores{i_exp}.basefit        = GLMType.fitname; 
    agg_scores{i_exp}.rawscores_note = 'rows are cells, columns [nullmodel, fit50, full-fit]';
    agg_scores{i_exp}.finalscore_note = '(objval(full)-objval(50)) / (objval(50)-objval(null))';
    for i_celltype = [1,2]
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_WN   = NaN(length(cellgroup),3);
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM = NaN(length(cellgroup),3);
        
        agg_scores{i_exp}.celltype{i_celltype}.finalscore_WN   = NaN(length(cellgroup),1);
        agg_scores{i_exp}.celltype{i_celltype}.finalscore_NSEM = NaN(length(cellgroup),1);
    end 
    
    % Assign Raw Scores
    for i_celltype = 1:2
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_WN(:,1) = objval.fit0.celltype{i_celltype}.scores_WN;
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_WN(:,2) = objval.fit50.celltype{i_celltype}.scores_WN;
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_WN(:,3) = objval.fitFULL.celltype{i_celltype}.scores_WN;
  
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM(:,1) = objval.fit0.celltype{i_celltype}.scores_NSEM;
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM(:,2) = objval.fit50.celltype{i_celltype}.scores_NSEM;
        agg_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM(:,3) = objval.fitFULL.celltype{i_celltype}.scores_NSEM;
    end
    
    for i_celltype = 1:2
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        
        for i_cell = 1:length(cellgroup)
            vec_WN = agg_scores{i_exp}.celltype{i_celltype}.rawscores_WN(i_cell,:);
            vec_NSEM = agg_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM(i_cell,:);
            
            if strcmp(shortname, 'halfratio_standardrk1')  
                if sum(isnan(vec_WN)) == 0
                    agg_scores{i_exp}.celltype{i_celltype}.finalscore_WN(i_cell) = ...
                        ( (vec_WN(3)-vec_WN(2)) / (vec_WN(2) - vec_WN(1)) );
                end
                if sum(isnan(vec_NSEM)) == 0
                    agg_scores{i_exp}.celltype{i_celltype}.finalscore_NSEM(i_cell) = ...
                        ( (vec_NSEM(3)-vec_NSEM(2)) / (vec_NSEM(2) - vec_NSEM(1)) );
                end
            else
                error('need to define new conv metric computation')
            end
            
            
        end
    end
end

aggregated_scores = agg_scores;
eval(sprintf('save %s/%s.mat aggregated_scores', Dirs.convergence, shortname));

end
%% PLOT
if runoptions.plot
 
eval(sprintf('load %s/%s.mat aggregated_scores', Dirs.convergence, shortname));


if strcmp(shortname, 'halfratio_standardrk1')   
    low_lim  = 0;
    high_lim = .2;
else
    error('need to register new conv metric for plotting') 
end
colors = {'r.','g.','b.','c.'};

for i_exp = 1:4
    clf
    colorstring = colors{i_exp};
    exp_nm = allcells{i_exp}.exp_nm;
    
    subplot(3,1,1);
    axis off; delta_y = .15;
    set(gca, 'fontsize', 12)
    c = 0; 
    if strcmp(shortname, 'halfratio_standardrk1')
        c=c+1; text(-.1, 1-delta_y*c,sprintf('Convergence Metric defined by objective value of fit set.  -- Named: %s', shortname),...
            'interpreter','none');
        c=c+1; text(-.1, 1-delta_y*c,'Metric Formula: (obval(fullfit)-objval(50pctfit)) / (objval(50pctfit) - objval(nullmodel))',...
            'interpreter','none');
    end
    c=c+1; text(-.1, 1-delta_y*c, sprintf('Model Fit: %s', GLMType.fitname),'interpreter','none' );
    c=c+1; text(-.1, 1-delta_y*c, sprintf('Plot Date %s',datestr(clock)),'interpreter','none');
    c=c+1; text(-.1, 1-delta_y*c, sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none' );
    
    for i_celltype = [1:2]
        
        if i_celltype == 1, subplot(3,2,[3 5]); end
        if i_celltype == 2, subplot(3,2,[4 6]); end

        MS_A = 12;
        MS_B = 30; 
        hold on; set(gca, 'fontsize', 12)
        axis square
        xlim([low_lim, high_lim]);
        ylim([low_lim, high_lim]);
        plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');


        scores1 = aggregated_scores{i_exp}.celltype{i_celltype}.finalscore_WN;
        scores2 = aggregated_scores{i_exp}.celltype{i_celltype}.finalscore_NSEM;
        title('Improvement 2nd half / Improvement 1st half');


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

        xlabel('White Noise')
        ylabel('NSEM')

    end    
    orient landscape
    eval(sprintf('print -dpdf %s/%s_%s.pdf',Dirs.convergence, shortname,exp_nm))
end

end



end        


function concat_fitmovie   = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , FitPars)
% AKHeitman 2014-04-14
% Concatenate the fit movie (different blocks)
% FitPars needs
%   .width
%   .height
%   .FitBlocks
%   .novelblocks
%   .fitframes

height       = FitPars.height;
width        = FitPars.width;
fitblocks    = FitPars.FitBlocks;
fitframes    = FitPars.fitframes;
novelblocks  = FitPars.NovelBlocks;

fitframesperblock = length(fitframes) ;
totalframes       = length(fitblocks) * ( fitframesperblock) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = fitblocks
        blkind = find(fitblocks == i_blk);
        framenums = ( (blkind -1)*fitframesperblock + 1 ) :  (blkind *fitframesperblock);  
        n_blkind = find(novelblocks == i_blk);
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{n_blkind}.matrix (:,:, fitframes);    
end

concat_fitmovie = concat_fullfitMovie;

end
  
function spikesconcat      = subR_concat_fitspikes_fromorganizedspikes(blockedspikes, FitPars)
% AKHeitman 2014-04-14
% Concatenate Spikes from different blocks to a single spike train
% blocekdspikes: needs
%   .t_sp_withinblock
%
% FitPars needs
%   .fittest_skipseconds
%   .tstim
%   .fitframes
%   .FitBlocks


t_start   = FitPars.fittest_skipseconds;
tstim     = FitPars.computedtstim;
fitframes = FitPars.fitframes;
FitBlocks = FitPars.FitBlocks;


T_SP = []; blk_count = 0;
dur = tstim * length(fitframes);
for k = FitBlocks
	blk_count = blk_count + 1;
	t_sp_full = blockedspikes.t_sp_withinblock{k} ; % unit of time: sec, 0 for the onset of the block
	t_sp      = t_sp_full(find(t_sp_full >  t_start));
	t_sp = t_sp - t_start;
	t_spcontext = t_sp + ( blk_count -1 )*dur;
	T_SP = [T_SP ; t_spcontext];
end
spikesconcat = T_SP;
end
        
        
function [lcif_stim, obj_val] = subR_findobj_lcifstim(lcif_external, pstar, GLMType, GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins)
% AKHEITMAN 2015-07-15
% subRoutine which will get optimized for inding input NL

% Part 1: Find Stim driven lcif (with input
% Part 3: Add in external lcif and find obj_val 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Find Stim driven lcif (with input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord



if GLMType.CONVEX
   
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    pstar = pstar';
    lcif_stim = pstar(paramind.X) *X_bin_shift;
end


if ~GLMType.CONVEX
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    pstar = pstar;
    
    
    frame_shifts = 0:1:(GLMPars.stimfilter.frames-1);
    if isfield(GLMPars.stimfilter,'frames_negative')
        frame_shifts = -(GLMPars.stimfilter.frames_negative):1:(GLMPars.stimfilter.frames-1)*1;
    end
    pixels  = size(X_frame,1);
    frames  = size(X_frame,2);
    
    
    % PARAMS TO GET INCORPORTED IN COVARIATE VEC 
    TimeFilter  = pstar(paramind.time1);
    SpaceFilter = pstar(paramind.space1);
    
    % FIND SPATIAL FILTER COVARIATE VEC (USING TIMEFILTER)
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    timeconvStim = zeros( pixels , (frames+length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % FIND TEMPORAL FILTER COVARIATE VEC (USING SPATIAL FILTER)
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    
    % STIMULUS ADDED TWICE / DIVIDE EACH COVARIATE VEC BY TWO
    stim_covariate_vec =  [.5*spatial_covariatevec; .5*temporal_covariatevec];    
    lcif_stim = (pstar(paramind.X))'  * stim_covariate_vec;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Add in external lcif and find obj_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    t_bin        = t_bin;
    home_sptimes = fitspikes.home';
    home_spbins  = ceil(home_sptimes / t_bin);
    home_spbins  = home_spbins(find(home_spbins < bins) );
    total_lcif = lcif_stim + lcif_external;
    cif        = exp(total_lcif);
    obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));
end

end

function [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins)
% AKHEITMAN 2015-07-15  extract non-stim components of lcif


% Construct PS_Basis
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
    

    if isfield(GLMType, 'special_arg') && isfield(GLMType.special_arg,'PS_Constrain')
            ps_basis_0 = ps_basis; clear ps_basis
            v        = sum(ps_basis_0,1);
            v        = v / norm(v) ;
            orthog_v = null(v);
            COB      = [v', orthog_v] ;
            ps_basis = (inv(COB) * ps_basis_0')' ;
    end

            
end

% Find PS_bin and MU_bin (pre-param multiply)
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


% Multiply by param for final log-cif values
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
pstar = pstar';
total_lcif = 0;
if isfield(paramind, 'MU')
    lcif_nonstim.components.mu   = pstar(paramind.MU)*MU_bin;
    total_lcif = total_lcif + lcif_nonstim.components.mu;
end
if isfield(paramind, 'PS')
    lcif_nonstim.components.ps   = pstar(paramind.PS)*PS_bin;
    total_lcif = total_lcif +  lcif_nonstim.components.ps;
end
lcif_nonstim.total = total_lcif;



end


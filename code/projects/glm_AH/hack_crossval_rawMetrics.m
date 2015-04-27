% GOOD CODE!! WILL FINALIZE SOON AKHEITMAN 2015-03-18
% local function (fraction of variance) 
% RESET
%{
clear; close all; clc
celltypes    = [2]
exps         = [4]
stimtypes    = [2]
cell_subset          = 'glmconv_4pct';
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear'
changes_cell{2}.type = 'input_pt_nonlinearity';
changes_cell{2}.name = 'piece_linear_aboutmean';
hack_crossval_rawMetrics(exps,stimtypes,celltypes,changes_cell,cell_subset)
%}
function hack_crossval_rawMetrics(exps,stimtypes,celltypes,changes_cell,cell_subset, runoptions)
% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE

if exist('runoptions','var')
    if isfield(runoptions,'debug') && runoptions.debug
        debug = true;
    end
end





% SET DIRECTORIES / USING GLM becaause it has the raster.
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType)


pars.vkspd_msec         =  50; 
pars.fracvar_smoothbins = 12;


% Victor Bin 
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP AND COMPUTE
    
%%% REORGANIZE INTO A MORE APPROACHABLE STRUCTURE
for i_exp = exps
    for i_stimtype = stimtypes
        for i_celltype = celltypes
            % Housekeeping
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            secondDir.exp_nm    = exp_nm;
            secondDir.stim_type = stimtype;
            secondDir.map_type  = 'mapPRJ';
            secondDir.fitname   = GLMType.fitname;
            glmfitdir   = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
            savedir = sprintf('%s/CrossVal_RawMetrics/', glmfitdir);
            if ~exist(savedir,'dir'), mkdir(savedir); end
            map_type= 'mapPRJ';
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end

            if strcmp(cell_subset,'all')
                candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP]
            elseif strcmp(cell_subset,'shortlist'); 
                [~,candidate_cells,~]  = cell_list(i_exp, 'shortlist'); 
                candidate_cells = cell2mat(candidate_cells) ; 
            elseif strcmp(cell_subset,'glmconv_4pct')
                conv_column = 2; 
                conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
                conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
                candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
            end
            cellgroup = intersect(candidate_cells, cellgroup);
            cellgroup = fliplr(cellgroup)

            if exist('debug','var') && debug
                cellgroup = cellgroup(1:2);
            end
            %%
            for i_cell = 1:length(cellgroup)
                tStart = tic;
                clear fittedGLM cell_savename
                cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
                
                
                eval(sprintf('load %s/%s.mat', glmfitdir, cell_savename))
                % PULL OUT COMPONENTS NECESSARY
                glm_rast = fittedGLM.xvalperformance.rasters.glm_sim;
                rec_rast = fittedGLM.xvalperformance.rasters.recorded;
                raw_bps  = fittedGLM.xvalperformance.logprob_glm_bpspike;
                t_bin    = fittedGLM.t_bin;
                
                %%
                
                l2error_norm = loc_L2error(rec_rast,glm_rast, t_bin, pars.fracvar_smoothbins);
                if exist('debug','var') && debug
                    glm_rast = glm_rast(1:2,:);
                    rec_rast = rec_rast(1:2,:);
                end
                timescale   = pars.vkspd_msec / 1000;
                [vspkd_mean vspkd_std] = loc_rastvspkd_comp(rec_rast,glm_rast,timescale,t_bin);
                
                crossval_rawmetrics.cid = cid;
                crossval_rawmetrics.cell_savename = cell_savename;
                crossval_rawmetrics.celltype = celltype;
                crossval_rawmetrics.exp_nm = exp_nm;
                crossval_rawmetrics.fit_objval = fittedGLM.rawfit.objective_val;
                crossval_rawmetrics.metric_pars = pars;
                crossval_rawmetrics.l2error_normvar  = l2error_norm;
                crossval_rawmetrics.vspkd            = vspkd_mean;
                crossval_rawmetrics.bps              = raw_bps;
                crossval_rawmetrics.note_l2  = sprintf('L2 error / temporal variance (L2(glm-rast)/var(rast)) (smoothed with 10msec guassian)');
                crossval_rawmetrics.note_bps = sprintf('bps = Bits Per Spike:=(logprob(rast|glm)-logprob(rast|tonic))/spikecount');
                crossval_rawmetrics.note_vspkd = sprintf('VictorSpike/Spikecount with cost of 1 for 50 msecs');
                crossval_rawmetrics.note_fitobjval = sprintf('Not CrossValidated.. just a reference fit point');
                
                eval(sprintf('save %s/%s.mat crossval_rawmetrics', savedir, cell_savename))
                duration = toc(tStart);
                display(sprintf('Finished: %s %s %s: Computation Time in Minutes', expname,cell_savename, stimtype,duration/60));
                
            end
           
            
         
        end

    end
end

end


function [l2_error_normed] = loc_L2error(rastRec,rastSim, t_bin, smoothbins)
% AKHeitman 2014-11-24
% rastRec is a binary binned matrix of spikes recorded
% rastSim is the corresponding simulation designed to immitate
% t_bin is the size of each bin in the raster in seconds
% smoothbins is the smoothing time scale in bins 

% Redone AKHeitman 2015-03-18

if (size(rastRec,1) ~= size(rastSim,1)) || (size(rastRec,2) ~= size(rastSim,2))
    error('rasters need to be of the same size!')
end
reps           = size(rastRec,1);
sigma_bin      = smoothbins;


convolve_index  = [-4*sigma_bin:4*sigma_bin];
convolve_vec    = normpdf(convolve_index,0,sigma_bin) / sum(normpdf(convolve_index,0,sigma_bin) );


ratebinRec      = sum(rastRec,1);
ratebinSim      = sum(rastSim,1);

ratesmoothRec   = conv(ratebinRec, convolve_vec) / reps;
ratesmoothSim   = conv(ratebinSim, convolve_vec) / reps;
ratesmoothRec   = ratesmoothRec( (4*sigma_bin+1):(end-4*sigma_bin));
ratesmoothSim   = ratesmoothSim( (4*sigma_bin+1):(end-4*sigma_bin));

signal_variance    = var(ratesmoothRec);

l2_error           = (1 / size(rastRec,2)) * sum ( (ratesmoothRec - ratesmoothSim).^2 );
l2_error_normed    =  l2_error / signal_variance; 
end
function [vspkd_mean vspkd_std] = loc_rastvspkd_comp(raster_A,raster_B,timescale,bindur)
%%% PURPOSE %%%
% How far is Raster B from Raster A. Normalized by spikecount of A.
% Scores scales from 0 to 2  ( .5 is good, 1 is bad, 2 is utter failure)
%%% NOTE %%%
% for really different rasters, we compute a failure of 2
% Metric only of use within the correct range of spikes

%%% INPUTS  %%%
% raster  matrix of zeros and ones
%         rows are repitions columns are time bins
% Time scale is in seconds
% bindur  is in seconds and is the duration of raster bin

%%% OUTPUTS %%%
% Crossval_Raw

%%% CALLS %%%
% loc_spkd (External Viktor Spike script)
%          Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%          Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.

%%%%%%%%%%%%%%%
% AKHeitman 2014-12-09 

% Version 1: Try to lower input param count
% Version 0: What worked a while ago in 2014
%%

% BASIC QUANTITIES
reps        = size(raster_A,1); 
bins        = size(raster_A,2); 
cost_param  = 1/timescale;
% LOOP THROUGH THE DIFFERENT TIME SCALES

    % COMPUTE VIKTOR SPIKE METRIC / LOOP THROUGH  
vksp_perspike = zeros(reps,1);
for i_rep = 1:reps
	spt_1        = bindur * find(raster_A(i_rep,:));
	spt_2        = bindur * find(raster_B(i_rep,:));
        
    % IF SPIKE COUNT IS WAY OFF IT SHOULD BE MARKED AS A FAILURE
	if length(spt_2) > 2*length(spt_1) || length(spt_2) < .5*length(spt_1)
        vksp_perspike(i_rep) = 2;
    else
        vksp_perspike(i_rep) = loc_spkd(spt_1, spt_2, cost_param) / (length(spt_1));
	end
end
    
vspkd_mean = mean(vksp_perspike);
vspkd_std  = std(vksp_perspike);
end



function d=loc_spkd(tli,tlj,cost)
%
% d=spkd(tli,tlj,cost) calculates the "spike time" distance
% (Victor & Purpura 1996) for a single cost
%
% tli: vector of spike times for first spike train
% tlj: vector of spike times for second spike train
% cost: cost per unit time to move a spike
%
%  Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%  Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
%
nspi=length(tli);
nspj=length(tlj);

if cost==0
   d=abs(nspi-nspj);
   return
elseif cost==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
%
%     INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
%
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if nspi & nspj
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+cost*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);

end
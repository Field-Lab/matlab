%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AKHeitman 2015-04-02

% Creates structure which dictates GLMType
% Loads cells 
% Loads stimuli / basic stimuli processing
% Loads spike trains / basic spike train processing
% Requires the organizedspikes structure with spike times relative
%    to start of each block of stimulus
% No direct GLM Paramater usage
% Feeds into glm_execute which is located in glm_core directory
% glm_execute along with glm_core 
%    which has no additional code dependencies, no loading of matfiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wrap_bookeeping Calls %
%  NSEM_BaseDirectories
%  GLM_Settings
%  GLM_fitname
%  NSEM_secondaryDirectories
%  loadmoviematfiles
%  StimulusParams

% Main Call %
%   glm_execute  

% Subroutines at bottom of function
%  subR_concat_fitspikes_fromorganizedspikes
%  subR_createraster
%  subR_concat_fitmovie_fromblockedcell
%  subR_visionSTA_to_xymviCoord
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
clear; clc
exps = 1; stimtypes = 1; celltypes = 2; 
cell_subset = 'debug'; postfilterNL.debug = true;
base_glmsettings = {}
glm_postfilterNL_wrap(exps,stimtypes,celltypes,cell_subset,base_glmsettings,postfilterNL)

%}
function glm_postfilterNL_wrap(exps,stimtypes,celltypes,cell_subset,base_glmsettings,postfilterNL)

% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
% Define structure which uniquely defines GLM to be used 
if exist('glm_settings', 'var')
    base_GLMType = GLM_settings('default',base_glmsettings);
else
    base_GLMType = GLM_settings('default');
end
base_GLMType.fitname    = GLM_fitname(base_GLMType); 


for i_exp = exps    
    for i_stimtype = stimtypes
        % Load master datarun, bookkeep
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load %s/%s/datarun_master.mat', BD.BlockedSpikes,exp_nm));
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        base_GLMType.fit_type = stimtype;
        
        % Load and process stimulus
        [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, base_GLMType.map_type);
        [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , stimtype, base_GLMType.cone_model,'fitmovie');
        [testmovie0]          = loadmoviematfile(exp_nm , stimtype, base_GLMType.cone_model,'testmovie');
        testmovie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
        base_GLMType.fitmoviefile  = origmatfile;
        if postfilterNL.debug
            StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:2);
        end
        fitmovie_concat       = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv); 
         
        % Directories  
        secondDir.exp_nm    = exp_nm; 
        secondDir.map_type  = base_GLMType.map_type; 
        secondDir.stim_type = stimtype;
        secondDir.fitname   = base_GLMType.fitname;
        Dirs.WN_STAdir          = NSEM_secondaryDirectories('WN_STA', secondDir); 
        Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
        Dirs.baseglm   = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        for i_celltype = celltypes            
            % Choose which subset of cells to run
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            if strcmp(cell_subset,'all')
                candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP]
            elseif strcmp(cell_subset,'shortlist') || strcmp(cell_subset, 'debug') 
                [~,candidate_cells,~]  = cell_list(i_exp, cell_subset); 
                candidate_cells = cell2mat(candidate_cells) ; 
            elseif strcmp(cell_subset,'glmconv_4pct')
                eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection));              
                conv_column = 2; 
                conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
                conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
                candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
            end
            cellgroup = intersect(candidate_cells, cellgroup)
        
            for i_cell = 1:length(cellgroup)
                cid = cellgroup(i_cell); 
                cell_savename = sprintf('%s_%d', celltype,cid);             
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                % Process spikes for glm_execute with proper subroutines
                fitspikes_concat.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
                
                % load fittedGLM
                eval(sprintf('load %s/%s.mat fittedGLM', Dirs.baseglm, cell_savename));
                glm_cellinfo = fittedGLM.cellinfo;
                pstar = fittedGLM.rawfit.opt_params;
                [lcif, obj_val] = decompose_covariates(pstar,fittedGLM.GLMType,fitspikes_concat,fitmovie_concat,inputstats,glm_cellinfo);  

                display(sprintf('For %s: %s: ObjVal is %d', expname, cell_savename, obj_val));
                % Loop here through various types of non-linearities
                t_bin = fittedGLM.t_bin;
                bins =  length(lcif.mu);
                home_sptimes = fitspikes_concat.home';
                home_spbins  = ceil(home_sptimes / t_bin);
                home_spbins  = home_spbins(find(home_spbins < bins) );
                optim_struct = optimset(...
                    'derivativecheck','off',...
                    'diagnostics','off',...  % 
                    'display','iter',...  %'iter-detailed',... 
                    'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
                    'GradObj','on',...
                    'largescale','on',...
                    'Hessian','on',...
                    'MaxIter',100,... % you may want to change this
                   'TolFun',10^(-(6)),...
                   'TolX',10^(-(9))   );
                if strcmp(postfilterNL.type, 'Null') 
                    p0  = [1 1 1]';   % should come back to [1 1 1]
                    COV = [lcif.stim; lcif.mu; lcif.ps];
                    [pstar fstar eflag output] = fminunc(@(p) subr_optimizationfunction...
                        (p,COV,home_spbins,t_bin),p0,optim_struct);
                    pstar
                elseif strcmp(postfilterNL.type, 'Hinge')
                    stim_pos = lcif.stim;
                    stim_pos(find(stim_pos<0)) = 0;
                    stim_neg = lcif.stim;
                    stim_neg(find(stim_neg>=0)) = 0;
                    
                    p0 = [1 1 1 1]';
                    COV = [stim_pos; stim_neg; lcif.mu; lcif.ps];
                    [pstar fstar eflag output] = fminunc(@(p) subr_optimizationfunction...
                        (p,COV,home_spbins,t_bin),p0,optim_struct);
                    display(sprintf('Pos Stim scale: %d, Neg Stim scale: %d', pstar(1),pstar(2)))
                    pstar
                end
            end            
        end
    end
end

end


function [f grad Hess log_cif]= subr_optimizationfunction(linear_params,covariates,spikebins,bin_duration)
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;
% Find Conditional Intensity and its log
lcif = p' * COV;
cif  = exp(lcif);
% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);
% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );
% Evaluate the hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');
% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
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
function raster_spiketimes = subR_createraster(blockedspikes, TestPars)
% AKHeitman 2014-04-14
% Make a raster which takes into account GLM processing
% blocekdspikes: needs
%   .t_sp_withinblock
%
% TestPars needs
%   .fittest_skipseconds
%   .TestBlocks

rasterblocks = TestPars.TestBlocks;
t_start      = TestPars.fittest_skipseconds;

raster_spiketimes = cell(length(rasterblocks),1);

for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = blockedspikes.t_sp_withinblock{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    % HACK NEEDED FOR 2013-10-10-0 and other long runs
    if isfield(TestPars, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (TestPars.test_skipENDseconds - TestPars.fittest_skipseconds - .1)));
    end
    
    raster_spiketimes{i_blk} = sptimes;
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
function [center,sd]       = subR_visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, masterdim, slvdim)
% AKHeitman  2013-12-08
% Grab x, y coordinates of STA center of the master
% Convert to coordinates of the enslaved dataset 
x_coord   = round( stafit_centercoord(1)* (slvdim.width  /masterdim.width)  );
y_coord   = slvdim.height - round( stafit_centercoord(2)* (slvdim.height /masterdim.height) );

center.x_coord = x_coord;
center.y_coord = y_coord;

sd.xdir = round( stafit_sd(1)* (slvdim.width   / masterdim.width)  );
sd.ydir = round( stafit_sd(2)* (slvdim.height  / masterdim.height)  );

end




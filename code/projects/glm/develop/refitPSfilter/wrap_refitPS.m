% 2015-07-02
% Specifically for optimized NLN model %
%{
clear; cell_subset = 'glmconv_4pct';
newPS = 'netinhibitory_domainconstrain_COB';
baseGLM.settings{1}.type = 'cone_model';
baseGLM.settings{1}.name = 'rieke_linear';
baseGLM.settings{2}.type= 'input_pt_nonlinearity';
baseGLM.settings{2}.name= 'piecelinear_fourpiece_eightlevels';
baseGLM.special_arg = 'Logistic_fixMU_noPS';
exps = [1 2 3 4];
stimtypes = [1 2];
celltypes = [1 2];
wrap_refitPS(baseGLM, newPS, exps, stimtypes, celltypes, cell_subset)
%}
function wrap_refitPS(baseGLM, newPS, exps, stimtypes, celltypes, cell_subset)
%%%%%%%%%
% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));

baseGLM.Type = GLM_settings('default',baseGLM.settings);
baseGLM.Type.fitname = GLM_fitname(baseGLM.Type); 
if isfield(baseGLM, 'special_arg')
    baseGLM.Type.fitname = sprintf('%s/%s', baseGLM.Type.fitname, baseGLM.special_arg);
end
   GLMType = baseGLM.Type;
GLMType.fitname   = sprintf('%s/refitPS_%s', GLMType.fitname, newPS);
GLMType.newPS  = newPS; 
GLMType.func_sname    = 'glm_wrap_constrainPS';
GLMType.fullmfilename = mfilename('fullpath'); 
display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));


% Run options, order cells for fitting
if exist('runoptions','var')
    if isfield(runoptions,'replace_existing')
        replace_existing  = true;
    end
    if isfield(runoptions,'reverseorder') 
        reverseorder  = true;
    end
end

for i_exp = exps    
    for i_stimtype = stimtypes
        % Load master datarun, bookkeep
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load %s/%s/datarun_master.mat', BD.BlockedSpikes,exp_nm));
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        GLMType.fit_type = stimtype;
        

        [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
        % Loading and Managing Directories
        secondDir.exp_nm        = exp_nm; 
        secondDir.map_type      = 'mapPRJ'; 
        secondDir.stim_type     = stimtype;
        secondDir.fitname       = baseGLM.Type.fitname;
        Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
        Dirs.baseglm_loaddir            = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        
        % Hack to get the correct save directory  
        BD_hack = BD;
        BD_hack.GLM_output_raw = BD.GLM_develop_output_raw;
        secondDir.fitname = GLMType.fitname;
        Dirs.fittedGLM_savedir  = NSEM_secondaryDirectories('savedir_GLMfit', secondDir,'',BD_hack)
        if ~exist(Dirs.fittedGLM_savedir), mkdir(Dirs.fittedGLM_savedir); end                  
        display(sprintf('Save Directory :  %s', Dirs.fittedGLM_savedir));
                
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
            cellgroup = intersect(candidate_cells, cellgroup);
            if exist('reverseorder','var') && reverseorder, cellgroup = fliplr(cellgroup); end           
            
            for i_cell = 1:length(cellgroup)
                cid = cellgroup(i_cell); 
                cell_savename = sprintf('%s_%d', celltype,cid);      
                eval(sprintf('load %s/%s fittedGLM', Dirs.baseglm_loaddir,cell_savename)) 
                fittedGLM_prePSrefit = fittedGLM; clear fittedGLM
                
                glm_cellinfo = fittedGLM_prePSrefit.cellinfo;
                glm_cellinfo.fitname = GLMType.fitname;
                glm_cellinfo.d_save = Dirs.fittedGLM_savedir;
                glm_cellinfo.t_bin  = fittedGLM_prePSrefit.t_bin;
                
                
                stimdrivenrate.fit  = fittedGLM_prePSrefit.stimtransform.cif_withNL_fit;
                stimdrivenrate.test = fittedGLM_prePSrefit.stimtransform.cif_withNL_test;
                
                % Load Blocked-Spikes from preprocessing
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                fitspikes_concat.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
                
                % Hack just load the premade raster
                testspikes_logicalraster = fittedGLM_prePSrefit.xvalperformance.rasters.recorded;
                
                if ~exist(sprintf('%s/%s.mat', Dirs.fittedGLM_savedir,cell_savename)) || (exist('replace_existing','var') && replace_existing)
                    [fittedGLM] = glm_execute_refitPS(newPS, GLMType,stimdrivenrate, fitspikes_concat, testspikes_logicalraster,glm_cellinfo);
                end
                
            end
        end
    end
end
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



%{
lcif_stim = log(fittedGLM.stimtransform.cif_withNL_fit);
eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));                    
% Process spikes for glm_execute with proper subroutines
fitspikes_concat.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);

GLMPars = GLMParams;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if COB
    ps_basis_0 = ps_basis; clear ps_basis
    v        = sum(ps_basis_0,1);
    v        = v / norm(v) ;
    orthog_v = null(v);
    COB      = [v', orthog_v] ;
    ps_basis = (inv(COB) * ps_basis_0')' ;
end

basis         = ps_basis';
PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
%}
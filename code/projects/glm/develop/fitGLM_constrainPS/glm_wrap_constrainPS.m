% STARTED 2015-06-29
% Will try to hold down the fitted PS filter

% Output: standardparams/ps_constrain.type

%{

clear
exps = [3];
stimtypes = [1];
celltypes = [1];
cell_subset = 'glmconv_4pct';
glm_settings{1}.type = 'cone_model';
glm_settings{1}.name = 'rieke_linear'
glm_settings{2}.type= 'input_pt_nonlinearity';
glm_settings{2}.name= 'piecelinear_fourpiece_eightlevels';
ps_constrain.type = 'PS_netinhibitory_domainconstrain_COB';
glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings, ps_constrain)


clear
exps = [1 3 2 4];
stimtypes = [2];
celltypes = [2 1];
cell_subset = 'glmconv_4pct';
glm_settings = {};
ps_constrain.type = 'PS_netinhibitory_domainconstrain_COB';
glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings, ps_constrain)

clear
exps = [1 3 2 4];
stimtypes = [1];
celltypes = [2 1];
cell_subset = 'all';
glm_settings = {};
ps_constrain.type = 'PS_netinhibitory_domainconstrain_COB';
glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings, ps_constrain)




clear
exps = [1 3 2 4];
stimtypes = [1];
celltypes = [2 1];
cell_subset = 'all';
glm_settings = {};
ps_constrain.type = 'PS_netinhibitory_domainconstrain_COB';
glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings, ps_constrain)



clear
exps = [1 2];
stimtypes = [2 1];
celltypes = [1 2];
cell_subset = 'all';
glm_settings = {};
%ps_constrain.type = 'PS_inhibitorydomainconstrain_post10msec';
ps_constrain.type = 'PS_netinhibitory_domainconstrain';
glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings, ps_constrain)

save('dbug_glmexecute_constrainPS', 'GLMType','fitspikes_concat','fitmovie_concat','testspikes_raster','testmovie','inputstats','glm_cellinfo')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function glm_wrap_constrainPS(exps,stimtypes,celltypes,cell_subset,glm_settings,ps_constrain,runoptions)

% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));

% Define structure which uniquely defines GLM to be used 
if exist('glm_settings', 'var')
    GLMType = GLM_settings('default',glm_settings);
else
    GLMType = GLM_settings('default');
end
baseGLMType = GLMType;
baseGLMType.fitname = GLM_fitname(GLMType); 

GLMType.fitname_preconstrainPS = GLM_fitname(GLMType); 
GLMType.fitname                = sprintf('%s/%s', GLMType.fitname_preconstrainPS,ps_constrain.type);
GLMType.PS_constrain   = ps_constrain; 

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
        
        % Load and process stimulus
        [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
        [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'fitmovie');
        [testmovie0]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
        testmovie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
        GLMType.fitmoviefile  = origmatfile;
        if GLMType.debug
            StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:2);
        end
        fitmovie_concat       = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv); 
         
        
        % Loading and Managing Directories
        secondDir.exp_nm        = exp_nm; 
        secondDir.map_type      = 'mapPRJ'; 
        secondDir.stim_type     = stimtype;
        secondDir.fitname       = baseGLMType.fitname;
        Dirs.WN_STAdir          = NSEM_secondaryDirectories('WN_STA', secondDir); 
        Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
        Dirs.baseglm            = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        
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
            cellgroup = intersect(candidate_cells, cellgroup)
            if exist('reverseorder','var') && reverseorder, cellgroup = fliplr(cellgroup); end           
            
            for i_cell = 1:length(cellgroup)
                cid = cellgroup(i_cell); 
                cell_savename = sprintf('%s_%d', celltype,cid);             
                if ~exist(sprintf('%s/%s.mat', Dirs.fittedGLM_savedir,cell_savename)) || (exist('replace_existing','var') && replace_existing)
                    % Create cell information structure
                    cell_savename
                    glm_cellinfo.cid            = cid;
                    glm_cellinfo.exp_nm         = exp_nm;
                    glm_cellinfo.celltype       = celltype;
                    glm_cellinfo.cell_savename  = cell_savename;
                    glm_cellinfo.fitname        = GLMType.fitname;
                    glm_cellinfo.d_save         = Dirs.fittedGLM_savedir;
                    glm_cellinfo.computedtstim  = StimulusPars.slv.computedtstim;
                    
                    % Add WN-STA and slave coordinates to glm_cellinfo
                    eval(sprintf('load %s/STAandROI_%s.mat STAandROI', Dirs.WN_STAdir, cell_savename));
                    master_idx         = find(datarun_master.cell_ids == cid);
                    stafit_centercoord = ( datarun_master.vision.sta_fits{master_idx}.mean );
                    stafit_sd          = ( datarun_master.vision.sta_fits{master_idx}.sd   );
                    slvdim.height      = StimulusPars.slv.height; 
                    slvdim.width       = StimulusPars.slv.width; 
                    [center_coord,sd]  = subR_visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
                    glm_cellinfo.WN_STA = STAandROI.STA;
                    glm_cellinfo.slave_centercoord = center_coord;
                    
                    % Load Blocked-Spikes from preprocessing
                    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                    
                    % Process spikes for glm_execute with proper subroutines
                    fitspikes_concat.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                    testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
                    tStart = tic;
                    % Load a previous initializer
                    
                    loadmatfile = sprintf('%s/%s.mat', Dirs.baseglm, cell_savename)
                    if exist(loadmatfile)
                        display('loading p_init')
                        eval(sprintf('load %s/%s.mat fittedGLM', Dirs.baseglm, cell_savename));
                        glm_cellinfo.p_init = fittedGLM.rawfit.opt_params;
                        if strcmp(baseGLMType.fitname,...
                                'fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL_stimnonlin_piecelinear_fourpiece_eightlevels/Man_DoubleOpt_standardparams')
                            glm_cellinfo.GLMPars = fittedGLM.GLMPars;
                        end
                        fittedGLM_old = fittedGLM; clear fittedGLM;
                    end
                   
                    tStart = tic;
                    [fittedGLM] = glm_execute_domainconstrainPS(ps_constrain.type,GLMType,...
                        fitspikes_concat,fitmovie_concat,testspikes_raster,testmovie,inputstats,glm_cellinfo);            
                    duration = toc(tStart);
                    display(sprintf('### runtime of %1.1e minutes ###', duration/60)); clear tStart duration tic
                               
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




function [testmovie, res] = res_spikes_raster(exp,celltype)
threshold = 0.1;
stimtype = 'NSEM';

datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8_shortlist/standardparams/NSEM_mapPRJ/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];

% Get all fits of that cell type
matfiles=dir([datapath exp_names(exp,:) celltype '*.mat']);
n_cells=length(matfiles);

% Load one fittedGLM to get the experiment info
load([datapath exp_names(exp,:) matfiles(1).name]);
GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;
% Load and process stimulus
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[movie, inputstats, ~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
testmovie             = movie{1}.matrix(:,:,StimulusPars.slv.testframes);
% Directories
secondDir.exp_nm    = exp_nm;
secondDir.map_type  = GLMType.map_type;
secondDir.stim_type = stimtype;
secondDir.fitname   = GLMType.fitname;
Dirs.fittedGLM_savedir  = NSEM_secondaryDirectories('savedir_GLMfit', secondDir);
Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);

% Loop through the cells and find the residual spikes for each cell.
res.spikes = zeros(n_cells, 10*length(testmovie));
res.cells = zeros(n_cells, 1);
res.centers = zeros(n_cells, 2);
for i_cell = 1:n_cells
    % Load fittedGLM
    load([datapath exp_names(exp,:) matfiles(i_cell).name]);
    % Load up neighbor spikes
    if GLMType.CouplingFilters
        n_couplings=length(fittedGLM.cellinfo.pairs);
        for i_pair=1:n_couplings
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir,  fittedGLM.cellinfo.pair_savename{i_pair}));
            neighborspikes{i_pair} = subR_createraster(organizedspikes.block, StimulusPars.slv);
        end
    else
        neighborspikes = 0;
    end
    % Calculate the rate from GLM
    rate = glm_rate_raster(fittedGLM, testmovie,inputstats,neighborspikes);
    spikes = fittedGLM.xvalperformance.rasters.recorded;
    res.spikes(i_cell,:) = sum(logical(spikes) & (rate < threshold));
    res.cells(i_cell) = fittedGLM.cellinfo.cid;
    res.centers(i_cell,1) = fittedGLM.cellinfo.slave_centercoord.y_coord;
    res.centers(i_cell,2) = fittedGLM.cellinfo.slave_centercoord.x_coord;
end

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

function concat_movie   = subR_concat_movie(movie)
% AKHeitman 2014-04-14
% Concatenate the fit movie (different blocks)
% FitPars needs
%   .width
%   .height
%   .FitBlocks
%   .novelblocks
%   .fitframes

height       = size(movie{1}.matrix, 2);
width        = size(movie{1}.matrix, 1);
fitblocks    = length(movie);
frames_per_block = size(movie{1}.matrix, 3);
concat_movie = uint8(zeros(width, height, fitblocks*frames_per_block )) ;

for i_blk = 1:fitblocks
        framenums = ((i_blk - 1)*frames_per_block+1):(i_blk*frames_per_block);
        concat_movie(:,:,framenums) = movie{i_blk}.matrix(:,:, :);    
end

end



function res = population_movie(stimtype, moviename, i_exp, i_celltype)

switch i_exp
    case 1; exp_nm = '2012-08-09-3';
    case 2; exp_nm = '2012-09-27-3';
    case 3; exp_nm = '2013-08-19-6';
    case 4; exp_nm = '2013-10-10-0';
end

% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
eval(sprintf('load %s/%s/datarun_master.mat', BD.BlockedSpikes,exp_nm));

% Load core directories and all eligible cells
%GLMType = fittedGLM.GLMType;
%exp_nm = fittedGLM.cellinfo.exp_nm;
GLMType.map_type = 'mapPRJ';

% Load and process test movie
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[testmovie0, ~,~]          = loadmoviematfile(exp_nm , stimtype, '8pix_Identity_8pix','testmovie');
testmovie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);

% Process spikes for glm_execute with proper subroutines
% Choose which subset of cells to run
if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end

n_cells = length(cellgroup);

testmovie_frames_per_block = length(testmovie);

res.spikes = zeros(n_cells, testmovie_frames_per_block);
res.cid = zeros(n_cells, 1);

for i_cell = 1:length(cellgroup)
    
    cid = cellgroup(i_cell);
    master_idx         = find(datarun_master.cell_ids == cid);
    cell_savename = sprintf('%s_%d', celltype,cid);
    
    % Load Blocked-Spikes from preprocessing
    secondDir.exp_nm    = exp_nm;
    secondDir.map_type  = GLMType.map_type;
    secondDir.stim_type = stimtype;
    Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
    
    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
    
    % Process spikes for glm_execute with proper subroutines
    testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
    
    
%     spikes_concat = [];
%     for i = 1:length(testspikes_raster.home)
%         trial_spikes = testspikes_raster.home{i};
%         spikes_concat = [spikes_concat; trial_spikes];
%     end
    spikes_frame = floor(testspikes_raster.home{1} * 240);
    for i_frame = 1:testmovie_frames_per_block
        res.spikes(i_cell, i_frame) = sum(spikes_frame == i_frame);
    end     
%     res.centers(i_cell,:) = flip(datarun_master.vision.sta_fits{master_idx}.mean);
    res.cid(i_cell) = cid;
end
res.spikes(res.spikes>1) = 1;
cell_spikes_movie(testmovie, res, moviename, datarun_master)

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


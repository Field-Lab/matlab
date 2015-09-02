function [spikes, movie, STA] = load_spikes_filter()

% load GLM to get info to load up movie and get STA
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
STA = flip(fittedGLM.cellinfo.WN_STA, 3);

% load spikes
load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/WN_mapPRJ/organizedspikes_ONPar_841.mat');

% Load core directories and all eligible cells
GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;
stimtype = fittedGLM.GLMType.fit_type;

% Load and process test movie
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[testmovie0, ~, ~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
movie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
clear testmovie0

% Process spikes for glm_execute with proper subroutines
spikes = subR_createraster(organizedspikes.block, StimulusPars.slv);

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
function [spikes, movie, STA] = load_spikes_filter(stimtype, exp_nm, cell)

% all inputs are strings
% stimtype is 'WN' or 'NSEM'
% exp_nm is '2012-08-09-3'
% cell is 'ONPar_841'

spike_directory = '/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/';

% load GLM to get info to load up movie and get white noise STA
load([spike_directory exp_nm '/WN_mapPRJ/STA/STAandROI_' cell '.mat'])
STA = flip(STAandROI.STA, 3);
STA = STA - mean(STA(:));
clear STAandROI

% Load and process test movie
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, 'mapPRJ');
[testmovie0, ~, ~]          = loadmoviematfile(exp_nm , stimtype, '8pix_Identity_8pix','testmovie');
movie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
movie = movie - mean(movie(:));
clear testmovie0

% Process spikes for glm_execute with proper subroutines
load([spike_directory exp_nm '/' stimtype '_mapPRJ/organizedspikes_' cell '.mat']);
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
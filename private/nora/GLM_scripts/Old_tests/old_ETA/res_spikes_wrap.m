function [testmovie, rate, res] = res_spikes_wrap(fittedGLM, stimtype, conv_blocks)

GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;

% Load and process stimulus
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[movie, inputstats, ~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'fitmovie');
inputstats.range = 255;
testmovie   = subR_concat_movie(movie((conv_blocks+1):end));

% Directories
secondDir.exp_nm    = exp_nm;
secondDir.map_type  = GLMType.map_type;
secondDir.stim_type = stimtype;
secondDir.fitname   = GLMType.fitname;
Dirs.fittedGLM_savedir  = NSEM_secondaryDirectories('savedir_GLMfit', secondDir);
Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);

% Create cell information structure
eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir,  fittedGLM.cell_savename));


% NBCoupling 2014-04-20
if GLMType.CouplingFilters
    n_couplings=length(fittedGLM.cellinfo.pairs); % number of cells to couple to
    % loading the neighboring spikes to neighborspikes.home
    for i_pair=1:n_couplings
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir,  fittedGLM.cellinfo.pair_savename{i_pair}));
        neighborspikes{i_pair} = subR_createraster(organizedspikes.block, StimulusPars.slv);
    end
else
    neighborspikes = 0;
end
% end NBCoupling

rate = glm_rate(fittedGLM,testmovie,inputstats,neighborspikes);
spikes = fittedGLM.xvalperformance.rasters.recorded;

for i_bin = 1:size(rate,2)
    rate_bin = rate(:,i_bin);
    spikes_bin = spikes(:,i_bin);
    spike_rate = rate_bin(logical(spikes_bin));
    res(i_bin) = sum(spike_rate > threshold);
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


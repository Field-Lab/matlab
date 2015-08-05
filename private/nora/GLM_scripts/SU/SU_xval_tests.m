
function lcif= SU_xval_tests(stimtype)

% load original fit and spikes
% load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8prefilter/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/WN_mapPRJ/organizedspikes_ONPar_841.mat');

% Load core directories and all eligible cells
GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;

% Load and process test movie
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[testmovie0, inputstats,~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
testmovie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);

% Process spikes for glm_execute with proper subroutines
testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
fittedGLM.GLMType.Subunits = 0;
fittedGLM.GLMType.timefilter = 'fit';
fittedGLM.GLMType.contrast = 0;

% lcif{1} = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats);
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8fit/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
% mess with the filters
% fittedGLM.GLMType.input_pt_nonlinearity = 1;
% fittedGLM.GLMType.input_pt_nonlinearity_type = 'exp';
% fittedGLM.GLMType.input_pt_nonlinearity = 0;
% fittedGLM.GLMType.Subunits = 1;
% fittedGLM.GLMType.Subunit_NL = 'exp';
% temp = [-0.1 0.5 0.5 0.5 -0.1];
% fittedGLM.SU_filter = temp'*temp;
lcif{1} = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats);
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

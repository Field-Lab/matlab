
function xvalperformance = glm_wrap_SU_iter(stimtype, iter)

% load original fit and spikes
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
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

% Use the second iteration to predict
fittedGLM.SU_filter = reshape(fittedGLM.rawfit.iter{iter}.SU, [3,3]);
paramind = fittedGLM.rawfit.paramind;
pstar = fittedGLM.rawfit.iter{iter}.nonSU;
basis_params  = fittedGLM.GLMPars.spikefilters.ps;
ps_basis      = prep_spikefilterbasisGP(basis_params,fittedGLM.t_bin);

% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
clear linearfilters
if isfield(paramind, 'MU')
    linearfilters.TonicDrive.Filter      = pstar(paramind.MU);
    linearfilters.TonicDrive.note        ='no convolution necessary, this term is part of the lcif for every bin';
end
if isfield(paramind, 'PS')
    rawfit.ps_basis = ps_basis;
    linearfilters.PostSpike.Filter     = ps_basis * pstar(paramind.PS);
    linearfilters.PostSpike.startbin   = 1;
    linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
    linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
end

% % SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
center_coord    = fittedGLM.cellinfo.slave_centercoord;
ROI_length      = fittedGLM.GLMPars.stimfilter.ROI_length;
stimsize.width  = size(testmovie,1);
stimsize.height = size(testmovie,2);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
% rawfit.ROIcoord = ROIcoord;
% clear stimsize center_coord;

ROI_length = fittedGLM.GLMPars.stimfilter.ROI_length;

if ~GLMType.CONVEX && (strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2'))
    if strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2')
        timefilter1  = pstar(paramind.time1);
        spacefilter1 = pstar(paramind.space1);
        stimfilter   = spacefilter1 * timefilter1';
        
        if strcmp(GLMType.stimfilter_mode, 'rk2')
            timefilter2  = pstar(paramind.time2);
            spacefilter2 = pstar(paramind.space2);
            stimfilter   = spacefilter1 * timefilter1' + spacefilter2 * timefilter2';
            
            % SVD to find Rank 1 and Rank 2 components
            [U,S,V]  = svd(reshape(stimfilter,[ROI_length^2,length(paramind.time1)] ));
            S = diag(S);
            xx = ( S(1)*U(:,1)*V(5,1) ) / norm( S(1)*U(:,1)*V(5,1) ) ;
            yy = S(1)*V(:,1) / norm( S(1)*V(:,1) );
            spacefilter1 = xx;
            timefilter1 = yy;
            
            xx = ( S(2)*U(:,2)*V(5,2) ) / norm( S(2)*U(:,2)*V(5,2) ) ;
            yy = S(2)*V(:,2) / norm( S(2)*V(:,2) );
            spacefilter2 = xx;
            timefilter2 = yy;
        end
        stimfilter = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.time1)]);
        linearfilters.Stimulus.Filter             = stimfilter;
        linearfilters.Stimulus.Filter_rank        = 1;
        linearfilters.Stimulus.time_rk1           = timefilter1;
        linearfilters.Stimulus.space_rk1          = reshape(spacefilter1,[ROI_length,ROI_length]);
        if strcmp(GLMType.stimfilter_mode, 'rk2')
            linearfilters.Stimulus.Filter_rank        = 2;
            linearfilters.Stimulus.time_rk2           = timefilter2;
            linearfilters.Stimulus.space_rk2          = reshape(spacefilter2,[ROI_length,ROI_length]);
        end
        linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
        linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
        linearfilters.Stimulus.frame_shifts       = 0:1:(fittedGLM.GLMPars.stimfilter.frames-1);
        linearfilters.Stimulus.bin_shifts         = 0:10:(fittedGLM.GLMPars.stimfilter.frames-1)*10;
        linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
        linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
        linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
        
        
    end
end

fittedGLM.linearfilters = linearfilters;
xvalperformance = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats,0);

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

%NBCoupling 2015-04-20
function paired_cells=subR_pick_neighbor_cells(mean, cell_ids, sta_fits)
    
     GLMPars = GLMParams;
     NumCells = length(cell_ids);
     distance=zeros(NumCells,1);
     
     % Calculate distance between RFs
     for i_pair=1:NumCells
         distance(i_pair)=norm(sta_fits{cell_ids(2,i_pair),1}.mean-mean);
         if distance(i_pair)==0
             distance(i_pair)=NaN;
         end
     end
     
     % Choose the closest cells
     [~,indices]=sort(distance);
     paired_cells=cell_ids(1,indices(1:GLMPars.spikefilters.cp.n_couplings));

end





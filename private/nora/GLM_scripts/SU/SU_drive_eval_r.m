iter = 2;

load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8prefilter/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
fittedGLM_orig = fittedGLM;
fittedGLM.SU_filter = reshape(fittedGLM_orig.rawfit.iter{2}.SU, [5,5]);
paramind = fittedGLM_orig.rawfit.paramind;
pstar = fittedGLM_orig.rawfit.iter{2}.nonSU;
GLMType = fittedGLM.GLMType;

%%

% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
clear linearfilters
if isfield(paramind, 'MU')
    linearfilters.TonicDrive.Filter      = pstar(paramind.MU);
    linearfilters.TonicDrive.note        ='no convolution necessary, this term is part of the lcif for every bin';
end
if isfield(paramind, 'PS')
    linearfilters.PostSpike.Filter     = fittedGLM.rawfit.ps_basis * pstar(paramind.PS);
    linearfilters.PostSpike.startbin   = 1;
    linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
    linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
end
% NBCoupling 05-28-14
if isfield(paramind, 'CP')
    rawfit.cp_basis = cp_basis;
    for j_pair=1:GLMPars.spikefilters.cp.n_couplings
        linearfilters.Coupling.Filter{j_pair}     = cp_basis * pstar(paramind.CP{j_pair});
    end
    linearfilters.Coupling.startbin   = 1;
    linearfilters.Coupling.note = 'Filter starts at "startbin" bins after the spikebin';
end
% end NBCoupling

% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
ROI_length      = fittedGLM.GLMPars.stimfilter.ROI_length;
ROIcoord = fittedGLM.rawfit.ROIcoord;
clear stimsize center_coord;

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
        linearfilters.Stimulus.frame_shifts       = 0:1:(GLMPars.stimfilter.frames-1);
        linearfilters.Stimulus.bin_shifts         = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
        linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
        linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
        linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';

    
    end  
end

lcif_SU = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats);
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/ONPar_841.mat')
lcif_noSU = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats);







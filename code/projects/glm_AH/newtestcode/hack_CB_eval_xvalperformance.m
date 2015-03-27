function [xvalperformance] = hack_CB_eval_xvalperformance(fittedGLM, SPars, organizedspikes,testmovie,inputstats)
%%% PURPOSE %%%
% Imitates eval_xvalperformance_NEW
% Specifically for Conductance Based Models
% CB Model fundamentally different because 2 distinct filters
% Evaluate Cross Validated Performance .. simulate fitted GLM

%%% INPUTS  %%%
% same as eval_xvalperformance_NEW

%%% OUTPUTS %%%
% xvalperformance
% same as eval_xvalperformance_NEW

% AKHEITMAN 2014-12-17
%% SAME AS EVAL_XVALPERFORMANCE_NEW
bpf = fittedGLM.bins_per_frame;
params.bindur = fittedGLM.t_bin;
params.bins = fittedGLM.bins_per_frame *length(SPars.testframes); 
params.evalblocks = SPars.TestBlocks;
params.trials = length(params.evalblocks);  
params.frames = length(SPars.testframes);
params.testdur_seconds = params.bindur * params.bins ;   

center_coord = fittedGLM.cellinfo.slave_centercoord;
teststim       = testmovie{1}.matrix;


frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord) ; 


logicalspike = zeros( length(params.evalblocks) , params.bins) ;         
for i_blk = 1 : length(params.evalblocks)
	blknum = params.evalblocks(i_blk);
	sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.fittest_skipseconds;
	sptimes = sptimes(find(sptimes > 0 ) );
    
    % NEW LINE TO ALLOW FOR SKIPPING BACK PART OF RASTER
    % HACK NEEDED FOR 2013-10-10-0 and other long runs
    % 2015-01-20
    if isfield(SPars, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (SPars.test_skipENDseconds - SPars.fittest_skipseconds - .1)));
    end
    
	for i_sp = 1:length(sptimes)
        spt = sptimes(i_sp);
        binnumber = ceil(spt / params.bindur );
        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
	end
end 
clear i_blk spt sptimes

GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same

[X_frame0 ] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, teststim,inputstats,center_coord) ;
X_frame     = X_frame0(:,SPars.testframes);

clear GLMType_fortest

GLMType = fittedGLM.GLMType;


 

MU = fittedGLM.linearfilters.TonicDrive.Filter;
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end


%% UNIQUE Conducatance Based Specific component
K_ex  = fittedGLM.linearfilters.Stimulus.Excitatory_Space * fittedGLM.linearfilters.Stimulus.Excitatory_Time';
K_in  = fittedGLM.linearfilters.Stimulus.Inhibitory_Space * fittedGLM.linearfilters.Stimulus.Inhibitory_Time';


KX_ex = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = prep_timeshift(X_frame(i_pixel,:),frame_shifts);
    tfilt = K_ex(i_pixel,:);
    KX_ex(i_pixel,:) = tfilt * X_frame_shift;
end
KX_ex0 = KX_ex;
lcif_kx_ex = sum(KX_ex,1);
lcif_kx_ex(find(lcif_kx_ex<=0)) = 0;

KX_in = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = prep_timeshift(X_frame(i_pixel,:),frame_shifts);
    tfilt = K_in(i_pixel,:);
    KX_in(i_pixel,:) = tfilt * X_frame_shift;
end
lcif_kx_in = sum(KX_in,1);
lcif_kx_in(find(lcif_kx_in>=0)) = 0;
%clf; plot(lcif_kx_ex(1:300)); hold on; plot(lcif_kx_in(1:300),'r'); hold off
lcif_kx_frame = lcif_kx_ex + lcif_kx_in;


%% SAME AS EVAL_XVALPERFORMANCE_NEW
lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
lcif_mu0 = MU * ones (1,params.bins);     
lcif_mu = repmat(lcif_mu0 , params.trials, 1);
lcif_kx = repmat(lcif_kx0 , params.trials, 1);    
clear sbpf;   
if GLMType.PostSpikeFilter
    lcif_ps = fastconvAH(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif_mu + lcif_kx + lcif_ps;
else
    lcif = lcif_mu + lcif_kx;
end
glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;
    
spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
model_null0      = spikerate_bin * ones(1, params.bins);
model_uop0       = (1/params.trials) * sum(logicalspike,1);
model_null       = repmat(model_null0, params.trials, 1);
model_uop        = repmat(model_uop0, params.trials, 1);
null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
uop_logprob      = sum(eval_rasterlogprob(logicalspike, model_uop, 'binary', 'conditioned'));
% Check computations are correct % 
%null_logprob    = sum(eval_rasterlogprob(logicalspike, model_null0, 'notbinary', 'unconditioned'));
%uop_logprob     = sum(eval_rasterlogprob(logicalspike, model_uop0, 'notbinary', 'unconditioned'));    
uop_bits             = uop_logprob - null_logprob;
uop_bits_perspike    = uop_bits / (sum(model_null0));
uop_bits_persecond   = uop_bits / params.testdur_seconds;
[raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
glm_logprob       = sum(raster_logprob_bin);
glm_bits          = glm_logprob - null_logprob;
glm_bits_perspike = glm_bits / (sum(model_null0));
glm_bits_perbin   = glm_bits / params.bins;
glm_bits_persecond   = glm_bits / params.testdur_seconds;


xvalperformance.logprob_null_raw     = null_logprob;
xvalperformance.logprob_uop_raw      =  uop_logprob;
xvalperformance.logprob_glm_raw      =  glm_logprob;
xvalperformance.logprob_uop_bpspike  =  uop_bits_perspike;
xvalperformance.logprob_glm_bpspike  =  glm_bits_perspike;
xvalperformance.logprob_uop_bpsec    =  uop_bits_persecond;
xvalperformance.logprob_glm_bpsec    =  glm_bits_persecond;
xvalperformance.glm_normedbits       =  glm_bits_persecond / uop_bits_persecond;


lcif_const  = lcif_kx0 + lcif_mu0;
logical_sim = zeros(params.trials, params.bins);


if GLMType.PostSpikeFilter
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : size(logicalspike,1)
        cif0         = exp(lcif_const);         
        cif_ps       = cif0;
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins- ps_bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif_ps(i));
                cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
    
   
    
else
    for i_trial = 1 : size(logicalspike,1)
        cif         = exp(lcif_const);         
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif(i));
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
end
    

xvalperformance.rasters.recorded = logicalspike;
xvalperformance.rasters.glm_sim  = logical_sim;
xvalperformance.rasters.bintime  = params.bindur;


end






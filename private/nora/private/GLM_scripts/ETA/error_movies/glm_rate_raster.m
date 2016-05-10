
function [glm_rateperbin] = glm_rate_raster(fittedGLM,testmovie,inputstats,testneighborspikes)
%%
bpf               = fittedGLM.bins_per_frame;
params.bindur     = fittedGLM.t_bin;
params.bins       = fittedGLM.bins_per_frame *size(testmovie,3); 
params.frames     = size(testmovie,3);
params.testdur_seconds = params.bindur * params.bins ;   
center_coord = fittedGLM.cellinfo.slave_centercoord;
teststim     = testmovie;
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

% Recorded spikes
logicalspike = fittedGLM.xvalperformance.rasters.recorded;
params.trials = size(logicalspike, 1);

% NBCoupling 2015-04-20
if fittedGLM.GLMType.CouplingFilters
    for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
        pairspike{pair} = zeros(params.trials,params.bins) ;
        for i_blk = 1 : params.trials
            spt = testneighborspikes{pair}{i_blk};
            binnumber = ceil(spt / params.bindur );
            pairspike{pair}( i_blk, binnumber )  =  pairspike{pair}( i_blk,binnumber ) + 1;
        end
        clear i_blk spt sptimes
    end
end

%%
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
inputstats.mu_avgIperpix = 64;
inputstats.range = 255;
[X_frame] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, teststim,inputstats,center_coord) ;
clear GLMType_fortest
GLMType = fittedGLM.GLMType;


    %% Set up CIF Components
MU = fittedGLM.linearfilters.TonicDrive.Filter;
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end
% NBCoupling 06-23-2014
if GLMType.CouplingFilters
    CP = fittedGLM.linearfilters.Coupling.Filter;
end
% NBcoupling
K  = fittedGLM.linearfilters.Stimulus.Filter;
K  = reshape(K, [ROI_pixels, length(frame_shifts)]);

KX = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = prep_timeshift(X_frame(i_pixel,:),frame_shifts);
    tfilt = K(i_pixel,:);
    KX(i_pixel,:) = tfilt * X_frame_shift;
end
lcif_kx_frame = sum(KX,1);

if isfield(GLMType, 'lcif_nonlinearity')
    lcif_kx_frame0 = lcif_kx_frame;
    
    if strcmp(GLMType.lcif_nonlinearity.type,'piece_linear_aboutmean')
        par = GLMType.lcif_nonlinearity.increment_to_decrement;
        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ;        
        pos_ind = find(lcif_kx_frame0>0);
        neg_ind = find(lcif_kx_frame0<=0);
        lcif_kx_frame = lcif_kx_frame0;
        lcif_kx_frame(pos_ind) = pos_mult * (lcif_kx_frame(pos_ind)); 
        lcif_kx_frame(neg_ind) = neg_mult * (lcif_kx_frame(neg_ind)); 
    elseif strcmp(GLMType.lcif_nonlinearity.type,'oddfunc_powerraise_aboutmean')
        par = GLMType.lcif_nonlinearity.scalar_raisedpower;       
        pos_ind = find(lcif_kx_frame0>0);
        neg_ind = find(lcif_kx_frame0<=0);
        lcif_kx_frame = lcif_kx_frame0;
        lcif_kx_frame(pos_ind) =  (     (lcif_kx_frame(pos_ind))  .*par );
        lcif_kx_frame(neg_ind) = -( (abs(lcif_kx_frame(neg_ind))) .*par );  
    end
end

%%
%display('binning the lcif components')
lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
lcif_mu0 = MU * ones (1,params.bins);     
lcif_mu = repmat(lcif_mu0 , params.trials, 1);
lcif_kx = repmat(lcif_kx0 , params.trials, 1);    
clear sbpf;  
lcif = lcif_mu + lcif_kx;
if GLMType.PostSpikeFilter
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif + lcif_ps; 
end
% NBCoupling 06-23-2014
if GLMType.CouplingFilters
    for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
        lcif_cp = fastconv(pairspike{pair} , [0; CP{pair}]', size(pairspike{pair},1), size(pairspike{pair},2) );
        lcif = lcif + lcif_cp;
    end
end
% end NBCoupling
glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;


end
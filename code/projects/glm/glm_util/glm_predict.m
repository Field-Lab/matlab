%% NB 2015-05-04

function [xvalperformance] = glm_predict(fittedGLM,testmovie,varargin)
%%
% INPUTS
% 
% fittedGLM structure
% testmovie should be in stim size x time (NO RGB!) and SAME STIXEL size as
%   the original fitmovie. It also gets rescaled, so don't worry about
%   scaling.
% OPTIONAL
% testspikes, which should be in cells, with each cell a repeat
%   if no testspikes, no bits per spike will be calculated
%   in SECONDS FROM BEGINNING OF REPEAT. This means you have to deal with
%   triggers, etc, BEFORE using in this function
% neighborspikes, if using coupling, same format as testspikes
% predict, set to 'false' if you dont want to make rasters, you just want
%   to calculate BPS for testspikes
% 

% Parse optional input 
p = inputParser;
p.addParamValue('testspikes', 0)
p.addParamValue('neighborspikes', 0)
p.addParamValue('predict', true)
p.addParamValue('trials', 20)
p.parse(varargin{:});
testspikes = p.Results.testspikes;
neighborspikes = p.Results.neighborspikes;
predict = p.Results.predict;
params.trials = p.Results.trials;
clear p


bpf               = fittedGLM.bins_per_frame;

if iscell(testspikes); params.trials     = length(testspikes); end % If there are testspikes, it will use that number of trials
params.bindur     = fittedGLM.t_bin;
params.bins       = fittedGLM.bins_per_frame *size(testmovie,3);
params.frames     = size(testmovie,3);
params.testdur_seconds = params.bindur * params.bins ;   
try
    center_coord = fittedGLM.center_coord;
catch
    center_coord = fittedGLM.cellinfo.slave_centercoord;
end
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

%%
if iscell(testspikes)
    logicalspike = zeros(params.trials,params.bins) ;
    for i_blk = 1 : params.trials
        spt = testspikes{i_blk};
        binnumber = ceil(spt / params.bindur );
        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
    end
    clear i_blk spt sptimes
    
end


% NBCoupling 2015-04-20
if fittedGLM.GLMType.CouplingFilters
    %try
        for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
            pairspike{pair} = zeros(params.trials,params.bins) ;
            for i_blk = 1 : params.trials
                if params.trials == 1
                    spt = neighborspikes{pair};
                else
                    spt = neighborspikes{pair}{i_blk};
                end
                binnumber = ceil(spt / params.bindur );
                pairspike{pair}( i_blk, binnumber )  =  pairspike{pair}( i_blk,binnumber ) + 1;
            end
            clear i_blk spt sptimes
        end
%     catch
%         warning('No neighbor spikes given, so coupling filters will be ignored.')
%         for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
%             pairspike{pair} = zeros(params.trials,params.bins) ;
%         end
%     end
end

%%
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
inputstats.mu_avgIperpix = mean(testmovie(:));
inputstats.range = max(testmovie(:))-min(testmovie(:));
[X_frame] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, testmovie,inputstats,center_coord) ;
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

if iscell(testspikes)
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
    xvalperformance.rasters.recorded = logicalspike;
    xvalperformance.glm_rateperbin  = params.bindur * glm_ratepersec;
    xvalperformance.glm_ratepersec = glm_ratepersec;
end

%%
if predict
    lcif_const  = lcif_kx0 + lcif_mu0;
    logical_sim = zeros(params.trials, params.bins);
    
    if GLMType.PostSpikeFilter && ~GLMType.CouplingFilters
        cif_psgain = exp(PS);
        ps_bins     = length(cif_psgain);
        for i_trial = 1 : size(logical_sim,1)
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
        
        
        % NBCoupling
    elseif GLMType.CouplingFilters && ~GLMType.PostSpikeFilter
        for i=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
            cif_cpgain{i} = exp(CP{i});
        end
        cp_bins     = length(cif_cpgain{1});
        for i_trial = 1 : size(logical_sim,1)
            cif0         = exp(lcif_const);
            cif_cp       = cif0;
            binary_simulation = zeros(1,params.bins);
            for i = 1 : params.bins- cp_bins;
                roll = rand(1);
                if roll >  exp(-params.bindur*cif_cp(i));
                    binary_simulation(i)= 1;
                end
                for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
                    if pairspike{pair}(i_trial,i)
                        cif_cp(i+1: i + cp_bins) =  cif_cp(i+1: i + cp_bins) .* (cif_cpgain{pair}');
                    end
                end
            end
            logical_sim(i_trial,:) = binary_simulation ;
            drive(i_trial, :) = params.bindur*cif_cp;
        end
        
    elseif GLMType.CouplingFilters && GLMType.PostSpikeFilter
        for i=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
            cif_cpgain{i} = exp(CP{i});
        end
        cif_psgain = exp(PS);
        ps_bins     = length(cif_psgain);
        cp_bins     = length(cif_cpgain{1});
        for i_trial = 1 : size(logical_sim,1)
            cif0         = exp(lcif_const);
            cif_ps_cp       = cif0;
            binary_simulation = zeros(1,params.bins);
            for i = 1 : params.bins- max(cp_bins, ps_bins);
                roll = rand(1);
                if roll >  exp(-params.bindur*cif_ps_cp(i));
                    cif_ps_cp(i+1: i + ps_bins) =  cif_ps_cp(i+1: i + ps_bins) .* (cif_psgain');
                    binary_simulation(i)= 1;
                end
                for pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
                    if pairspike{pair}(i_trial,i)
                        cif_ps_cp(i+1: i + cp_bins) =  cif_ps_cp(i+1: i + cp_bins) .* (cif_cpgain{pair}');
                    end
                end
            end
            logical_sim(i_trial,:) = binary_simulation ;
            drive(i_trial, :) = params.bindur*cif_ps_cp;
        end
        % end NBCoupling
        
        
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
    xvalperformance.rasters.glm_sim  = logical_sim;
end

xvalperformance.rasters.bintime  = params.bindur;


end






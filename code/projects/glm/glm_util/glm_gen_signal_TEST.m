%% NB 2015-05-04

function [lcif] = glm_gen_signal_TEST(fittedGLM,testmovie,varargin)
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

%display('binning the lcif components')
lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
lcif_mu0 = MU * ones (1,params.bins);
lcif = lcif_kx0+lcif_mu0;

end






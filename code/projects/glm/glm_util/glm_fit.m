%% NB 2015-05-01
% This takes in the stimulus, the spikes, and the location of the cell to
% fit a GLM. The GLM architecture and settings can be changed in
% glm_parameters.m
%
% PATHS NEEDED
% Vision.jar, such as javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar')
% The lab codebase, addpath(genpath('Repo location /matlab/code/lab'))
% The glm code folder, addpath(genpath('Repo location /matlab/code/projects/glm))

% TROUBLESHOOTING
% If it isn't working, use STA_Test(fitspikes, fitmovie, center_coord) to
% make sure your timing and indexing is right.
% This also only works for greyscale movies! It is not set up for color!


function [fittedGLM] = glm_fit(fitspikes, fitmovie, center, varargin)

% INPUTS
%

% REQUIRED

%   fitspikes: the spike times of the neuron

%   fitmovie: the movie frame by frame. You should
%   have a frame for every 1/120 seconds, so if the interval was two, your
%   fitmovie should have 2 of each frame
%   OR the xml specification, like RGB-8-1-0.48-11111-32x32

%   center_coord: the center of the RF

% OPTIONAL

%   WN_STA: optional, To do fixedSP_rk1, you need to input the STA in the same
%   dimensions as the fitting stimulus

%   neighborspikes: optional, a cell array, where each cell has the spike times of
%   the neighbor cells
%
%   monitor_refresh: usually should be 120Hz. This is NOT the interval!
%   Just the monitor speed!

% Parse optional input 
p = inputParser;
p.addParamValue('WN_STA', 0)
p.addParamValue('neighborspikes', 0)
p.addParameter('monitor_refresh', 120)
p.parse(varargin{:});
WN_STA = p.Results.WN_STA;
neighborspikes = p.Results.neighborspikes;
monitor_refresh = p.Results.monitor_refresh;
clear p

%% Setup Covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[GLMT, GLMP] = glm_parameters; % GLMType and GLMParams defined here

fittedGLM.GLMPars = GLMP;
fittedGLM.GLMType = GLMT;
GLMPars = GLMP;
GLMType = GLMT;
clear GLMP GLMT

center_coord.x_coord = center(2);
center_coord.y_coord = center(1);
fittedGLM.center_coord = center_coord;

if isfield(fittedGLM.GLMType, 'specialchange') && fittedGLM.GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end

if isfield(fittedGLM.GLMType, 'debug') && fittedGLM.GLMType.debug
    GLMPars.optimization.tolfun = 1; 
end

% Timing
frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = (1/monitor_refresh) / GLMPars.bins_per_frame;
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;

% Make Coupling and Post Spike Filter Bases
bin_size      = t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
clear bin_size basis_params

% Convolve Spike Times with appropriate basis
home_spbins  = ceil(fitspikes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
% NBCoupling 05-28-14
if GLMType.CouplingFilters;
    basis = cp_basis';
    for j_pair=1:GLMPars.spikefilters.cp.n_couplings
        neighbor_sptimes = neighborspikes{j_pair}';
        neighbor_spbins  = ceil(neighbor_sptimes / t_bin);
        neighbor_spbins = neighbor_spbins(find(neighbor_spbins < bins) );
        CP_bin{j_pair}=prep_convolvespikes_basis(neighbor_spbins,basis,bins);
    end
end
% end NBCoupling

if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end

% PREPARE PARAMETERS
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
p_init     =  0.1*ones(paramind.paramcount,1);  

% ORGANIZE STIMULUS COVARIATES
inputstats.mu_avgIperpix = double(mean(fitmovie(:)));
inputstats.range = double(range(fitmovie(:)));
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
fittedGLM.inputstats = inputstats;
%
if ~GLMType.CONVEX
    GLMPars.optimization.tolfun = GLMPars.optimization.tolfun - 1;
    fittedGLM.GLMPars = GLMPars;
end
%}
% INITIALIZE OPTIMIZATION STRUCTURE FOR MATLAB FMIN SOLVERS
if GLMType.CONVEX
    optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   );
end
if ~GLMType.CONVEX
    optim_struct = optimset(...
    'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'Hessian','on',...
   'largescale','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   );
end


%% Run through optimization .. get out pstart, fstar, eflag, output
% CONVEXT OPTIMIZATION
if GLMType.CONVEX
    glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(paramind, 'MU')
        glm_covariate_vec( paramind.MU , : ) = MU_bin;
    end
    if isfield(paramind, 'X')
        glm_covariate_vec( paramind.X , : ) = X_bin_shift;
    end
    if isfield(paramind, 'PS')
        glm_covariate_vec( paramind.PS , : ) = PS_bin;
    end
    % NBCoupling 05-28-14
    if isfield(paramind, 'CP')
        for j_pair=1:GLMPars.spikefilters.cp.n_couplings
            glm_covariate_vec( paramind.CP{j_pair} , : ) = CP_bin{j_pair};
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(GLMType, 'postfilter_nonlinearity') || ~GLMType.postfilter_nonlinearity
        [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),p_init,optim_struct);
    end
    if isfield(GLMType, 'postfilter_nonlinearity') && GLMType.postfilter_nonlinearity
        [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction_withNL...
            (p,glm_covariate_vec,home_spbins,t_bin,nonlinearity),p_init,optim_struct);
        % [f grad Hess log_cif COV_NL]=glm_convex_optimizationfunction_withNL(pstar,glm_covariate_vec,home_spbins,t_bin,nonlinearity);        
    end
end

% NONCONVEX OPTMIZATION
if ~GLMType.CONVEX
    bpf               = GLMPars.bins_per_frame;
    frame_shifts      = 0:1:(GLMPars.stimfilter.frames-1);
    % dnote part that is convex
    convex_cov = NaN(paramind.convParams , bins );  % convex covariates
    if isfield(paramind, 'MU')
        convex_cov( paramind.MU , : ) = MU_bin;
    end
    if isfield(paramind, 'PS')
        convex_cov( paramind.PS , : ) = PS_bin;
    end
    % NBCoupling
    if isfield(paramind, 'CP')
        for j_pair=1:GLMPars.spikefilters.cp.n_couplings
            convex_cov( paramind.CP{j_pair} , : ) = CP_bin{j_pair};
        end
    end
    filtertype = GLMType.stimfilter_mode;
    
    [pstar fstar eflag output] = fminunc(@(p) glm_nonconvex_optimizationfunction...
            (p,filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins,t_bin),p_init,optim_struct);

end
fittedGLM.fminunc_output = output;

%% Unpack the output into filters

rawfit.opt_params        = pstar;
rawfit.paramind          = paramind;
rawfit.objective_val     = fstar;

% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
clear linearfilters
linearfilters.note = 'These terms, convolved with the covariates, make the log of the conditional intensity function';
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
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
rawfit.ROIcoord = ROIcoord;
clear stimsize center_coord;
if GLMType.CONVEX
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        [STA_sp,STA_time]= spatialfilterfromSTA(WN_STA,ROIcoord.xvals,ROIcoord.yvals);
        timefilter           = pstar(paramind.X);
        stimfilter           = STA_sp * (timefilter');
        stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.X)]);
        rawfit.spatialfilter = STA_sp;
        linearfilters.Stimulus.Filter             = stimfilter;
        linearfilters.Stimulus.Filter_rank        = 1;
        linearfilters.Stimulus.space_rk1          = reshape(STA_sp, [ROI_length,ROI_length]);
        linearfilters.Stimulus.time_rk1           = pstar(paramind.X);
        %linearfilters.Stimulus.WN_note            = 'use WN STA as a reference to compare to fitted filters'
        %linearfilters.Stimulus.WN_STA             = WN_STA;
        %linearfilters.Stimulus.WN_STA_space_rk1   = reshape(STA_sp, [ROI_length,ROI_length]);
        %linearfilters.Stimulus.WN_STA_time_rk1    = STA_time;
        linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
        linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
        linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
        linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];
        
        if isfield(GLMPars.stimfilter,'frames_negative')
            linearfilters.Stimulus.frame_shifts       = [(-GLMPars.stimfilter.frames_negative):1:(GLMPars.stimfilter.frames-1)];
            linearfilters.Stimulus.bin_shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
        end
        linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
        linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
        linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
    end
end 
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
        linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
        linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];
        linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
        linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
        linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
    end  
end

fittedGLM.rawfit               = rawfit;
fittedGLM.linearfilters = linearfilters;
fittedGLM.note = 'in theory, linearfilters and t_bin/ binsperframe is sufficient for xval and simulation'; 
fittedGLM.fit_time = datestr(clock);
fittedGLM.writingcode = mfilename('fullpath');

end

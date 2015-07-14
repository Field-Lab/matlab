%% AKHeitman 2014-04-27
% This will be admittedly long and ugly.  But more self contained.
% GLMPars = GLMParams  or GLMPars = GLMParams(GLMType.specialchange_name)
% Sensitive to naming changes in GLMParams.
% Only saves stuff (and calls directories) if we are in a troubleshooting mode
% Heavily GLMType dependent computations will be carried out here
% Outsourcable computations will be made into their own functions
% troubleshoot optional
% need troublshoot.doit (true or false), 
% troubleshoot.plotdir,
% troubleshoot.name



% CALLS which use GLMType:
%  prep_paramindGP
%  prep_stimcelldependentGPXV


function [fittedGLM] = glm_fitcrossval_execute(GLMType,spikes_raster,movie,inputstats,glm_cellinfo,neighborspikes)

%% Setup Covariates
fittedGLM.cell_savename = glm_cellinfo.cell_savename;
fittedGLM.d_save        = glm_cellinfo.d_save;
fittedGLM.cellinfo      = glm_cellinfo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars           = GLMParams;
if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;
if isfield(GLMType, 'debug') && GLMType.debug
    GLMPars.optimization.tolfun = 1; 
end

% Timing
reps = length(spikes_raster.home);
frames = size(movie,3);
bins_perrep   = frames * GLMPars.bins_per_frame;
bins = bins_perrep * reps;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
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
    % Put in PCA for coupling here
    % load('CP_basis.mat');
    % cp_basis = waveform; 
end
clear bin_size basis_params

% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
t_bin        = t_bin;

home_spbins.perrep     = cell(reps,1);
home_spbins.continuous = []; 
for i_rep = 1:reps
    sptimes = spikes_raster.home{i_rep}';
    spbins  = ceil(sptimes / t_bin);
    spbins  = spbins(find(spbins < bins_perrep) );
    
    home_spbins.perrep{i_rep} = spbins;
    offset = (i_rep-1) * bins_perrep;
    offset_spbins = spbins + offset;
    home_spbins.continuous = [home_spbins.continuous,offset_spbins];    
end


if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end

% PREPARE PARAMETERS
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
%p_init     =  zeros(paramind.paramcount,1);  
p_init     = .01* ones(paramind.paramcount,1);


% ORGANIZE STIMULUS COVARIATES
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
if isfield(paramind, 'X')
    [X_frame_0,X_bin_0]    = prep_stimcelldependentGPXV(GLMType, GLMPars, movie, inputstats, center_coord, WN_STA);
    X_frame = repmat(X_frame_0, [1,reps]);
    X_bin   = repmat(X_bin_0, [1,reps]);
    
end



clear WN_STA center_coord

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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(paramind, 'MU')
        glm_covariate_vec( paramind.MU , : ) = MU_bin;
    end
    if isfield(paramind, 'X')
        X_bin_shift = prep_timeshift(X_bin,shifts);
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
        [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins.continuous,t_bin),p_init,optim_struct);
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
            (p,filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins.continuous,t_bin),p_init,optim_struct);

end
fittedGLM.fminunc_output = output;

%% Unpack the output into filters
testmovie = movie;
testspikes_raster = spikes_raster;
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
center_coord    = glm_cellinfo.slave_centercoord;
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(movie,1);
stimsize.height = size(movie,2);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
rawfit.ROIcoord = ROIcoord;
clear stimsize center_coord;
WN_STA           = double(glm_cellinfo.WN_STA); 
[STA_sp,STA_time]= spatialfilterfromSTA(WN_STA,ROIcoord.xvals,ROIcoord.yvals);
if GLMType.CONVEX
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        
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
    
    % Hacked no stim mode AH 2015-07-08
    if strcmp(GLMType.stimfilter_mode, 'nostim')
        
        timefilter           = zeros(GLMPars.stimfilter.frames,1);
        stimfilter           = STA_sp * (timefilter');
        stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,GLMPars.stimfilter.frames]);
        rawfit.spatialfilter = STA_sp;
        linearfilters.Stimulus.Filter             = stimfilter;
        linearfilters.Stimulus.Filter_rank        = 1;
        linearfilters.Stimulus.space_rk1          = reshape(STA_sp, [ROI_length,ROI_length]);
        linearfilters.Stimulus.time_rk1           = zeros(GLMPars.stimfilter.frames,1);
        %linearfilters.Stimulus.WN_note           = 'use WN STA as a reference to compare to fitted filters'
        %linearfilters.Stimulus.WN_STA             = WN_STA;
        %linearfilters.Stimulus.WN_STA_space_rk1   = reshape(STA_sp, [ROI_length,ROI_length]);
        %linearfilters.Stimulus.WN_STA_time_rk1    = STA_time;
        linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
        linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
        linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
        linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];

        linearfilters.Stimulus.note1              = 'Hack fill in so rest of code works.  Just put a zero stimulus fiter';
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

%% Evaluate cross-validated fits,  Print and Save
[xvalperformance] = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats,neighborspikes.test);
fittedGLM.xvalperformance  = xvalperformance; 
eval(sprintf('save %s/%s.mat fittedGLM',glm_cellinfo.d_save,glm_cellinfo.cell_savename));
printname = sprintf('%s/DiagPlots_%s',glm_cellinfo.d_save,fittedGLM.cellinfo.cell_savename);
printglmfit(fittedGLM,printname)

end

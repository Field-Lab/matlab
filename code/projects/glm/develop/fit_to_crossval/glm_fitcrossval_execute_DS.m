% AKHEITMAN 2015-07-10
% Exploratory code to see what happens if we fit with the crossval set
% Way to say that our fitting of linear filter isn't the problem
% CALLS which use GLMType:
%  prep_paramindGP
%  prep_stimcelldependentGPXV


function [fittedGLM] = glm_fitcrossval_execute_DS(crossval,GLMType,spikes_raster,movie,inputstats,glm_cellinfo,neighborspikes)

if strcmp(crossval.name, 'fit_crossval_DS')
    fitspikes_raster = spikes_raster;
    testspikes_raster = spikes_raster;
    fitmovie = movie;
    testmovie = movie;
end

if strcmp(crossval.name,'fit_crossval_oddeven_DS') || strcmp(crossval.name,'fit_crossval_oddeven_3DS')
    reps = floor(.5*length(spikes_raster.home));
    
    fitspikes_raster.home  = cell(reps,1);
    testspikes_raster.home = cell(reps,1); 
    for i_rep = 1:reps
        fitspikes_raster.home{i_rep}  = spikes_raster.home{2*i_rep-1};
        testspikes_raster.home{i_rep} = spikes_raster.home{2*i_rep};
    end
    
    fitmovie = movie;
    testmovie = movie;
end

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
reps = length(fitspikes_raster.home);
frames = size(fitmovie,3);
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
    sptimes = fitspikes_raster.home{i_rep}';
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
p_init     = .01* ones(paramind.paramcount,1);


% ORGANIZE STIMULUS COVARIATES
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
% DOWN SAMPLING

if strcmp(crossval.name,'fit_crossval_oddeven_DS')
    ROI_length0 = GLMPars.stimfilter.ROI_length;
    ROI_length  = floor(ROI_length0/2);
    [X_frame_0,X_bin_0]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
    
    xf = reshape(X_frame_0, [ROI_length0 ROI_length0 frames]);
    xb = reshape(X_bin_0, [ROI_length0 ROI_length0 bins_perrep]);
    if 2*ROI_length < ROI_length0
        xf = xf(2:end,2:end,:);
        xb = xb(2:end,2:end,:);
    end
    
    xf_ds = NaN([ROI_length ROI_length frames]);
    xb_ds = NaN([ROI_length ROI_length bins_perrep]);
    
    for i_row = 1:ROI_length
        for i_col = 1:ROI_length
            rows = 2*i_row + [-1 0];
            cols = 2*i_col + [-1 0];
            clump = xf([rows],[cols],:);
            new_vec = squeeze( mean(mean(clump,1),2) ); 
            xf_ds(i_row,i_col,:) = new_vec;
        end
    end
    for i_row = 1:ROI_length
        for i_col = 1:ROI_length
            rows = 2*i_row + [-1 0];
            cols = 2*i_col + [-1 0];
            clump = xb([rows],[cols],:);
            new_vec = squeeze( mean(mean(clump,1),2) ); 
            xb_ds(i_row,i_col,:) = new_vec;
        end
    end
    
    X_frameDS_0 = reshape(xf_ds, [ROI_length^2 frames]);
    X_binDS_0 = reshape(xb_ds, [ROI_length^2 bins_perrep]);
    
    X_frame = repmat(X_frameDS_0, [1,reps]);
    X_bin   = repmat(X_binDS_0, [1,reps]);  
    
    paramind_0 = paramind; clear paramind
    paramind.MU = 1;
    paramind.X  = 1 + [1:(ROI_length^2 + GLMPars.stimfilter.frames)];
    paramind.space1 = 1 + [1:ROI_length^2];
    paramind.time1  = 1 + ROI_length^2 + [1:GLMPars.stimfilter.frames];
    paramind.convParams = 1;
    paramind.convParams_ind = 1;
    paramind.paramcount = 1 + ROI_length^2 + GLMPars.stimfilter.frames;
    p_init     = .01* ones(paramind.paramcount,1);
end


if strcmp(crossval.name,'fit_crossval_oddeven_3DS')
    ROI_length0 = GLMPars.stimfilter.ROI_length;
    ROI_length  = floor(ROI_length0/3);
    [X_frame_0,X_bin_0]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
    
    xf = reshape(X_frame_0, [ROI_length0 ROI_length0 frames]);
    xb = reshape(X_bin_0, [ROI_length0 ROI_length0 bins_perrep]);
    if 3*ROI_length == (ROI_length0-1)
        xf = xf(2:end,2:end,:);
        xb = xb(2:end,2:end,:);
    end
    
    xf_ds = NaN([ROI_length ROI_length frames]);
    xb_ds = NaN([ROI_length ROI_length bins_perrep]);
    
    for i_row = 1:ROI_length
        for i_col = 1:ROI_length
            rows = 3*i_row + [-2 -1 0];
            cols = 3*i_col + [-2 -1 0];
            clump = xf([rows],[cols],:);
            new_vec = squeeze( mean(mean(clump,1),2) ); 
            xf_ds(i_row,i_col,:) = new_vec;
        end
    end
    for i_row = 1:ROI_length
        for i_col = 1:ROI_length
            rows = 3*i_row + [-2 -1 0];
            cols = 3*i_col + [-2 -1 0];
            clump = xb([rows],[cols],:);
            new_vec = squeeze( mean(mean(clump,1),2) ); 
            xb_ds(i_row,i_col,:) = new_vec;
        end
    end
    
    X_frameDS_0 = reshape(xf_ds, [ROI_length^2 frames]);
    X_binDS_0 = reshape(xb_ds, [ROI_length^2 bins_perrep]);
    
    X_frame = repmat(X_frameDS_0, [1,reps]);
    X_bin   = repmat(X_binDS_0, [1,reps]);  
    
    paramind_0 = paramind; clear paramind
    paramind.MU = 1;
    paramind.X  = 1 + [1:(ROI_length^2 + GLMPars.stimfilter.frames)];
    paramind.space1 = 1 + [1:ROI_length^2];
    paramind.time1  = 1 + ROI_length^2 + [1:GLMPars.stimfilter.frames];
    paramind.convParams = 1;
    paramind.convParams_ind = 1;
    paramind.paramcount = 1 + ROI_length^2 + GLMPars.stimfilter.frames;
    p_init     = .01* ones(paramind.paramcount,1);
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
   'MaxIter',100,... % you may want to change this
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
   'MaxIter',100,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   );
end


%% Run through optimization .. get out pstart, fstar, eflag, output

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
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
rawfit.ROIcoord = ROIcoord;

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
[xvalperformance] = hack_xvalperformance_DS(X_frameDS_0,fittedGLM,testspikes_raster,testmovie,inputstats,neighborspikes.test);
fittedGLM.xvalperformance  = xvalperformance; 
eval(sprintf('save %s/%s.mat fittedGLM',glm_cellinfo.d_save,glm_cellinfo.cell_savename));
printname = sprintf('%s/DiagPlots_%s',glm_cellinfo.d_save,fittedGLM.cellinfo.cell_savename);
printglmfit(fittedGLM,printname)

end

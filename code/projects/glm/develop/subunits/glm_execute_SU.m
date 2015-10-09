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


function [fittedGLM] = glm_execute_SU(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes)


%% NON STIMULUS STUFF, can leave alone, will be added in linearly 
% setup
fittedGLM.cell_savename = glm_cellinfo.cell_savename;
fittedGLM.d_save        = glm_cellinfo.d_save;
fittedGLM.cellinfo      = glm_cellinfo;
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
frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;

% Make bases for the other inputs
bin_size      = t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if isfield(GLMType,'Contrast') && GLMType.Contrast
    basis_params  = GLMPars.spikefilters.C;
    C_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if isfield(GLMType, 'Saccades')
   basis_params = GLMPars.saccadefilter;
   sa_basis = prep_spikefilterbasisGP(basis_params, bin_size);
end
clear bin_size basis_params

% Convolve with Bases to get covariates
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
center_coord       = glm_cellinfo.slave_centercoord;
if isfield(GLMType,'Contrast') && GLMType.Contrast
    stimsize.width  = size(fitmovie,1);
    stimsize.height = size(fitmovie,2);
    ROIcoord        = ROI_coord(GLMPars.spikefilters.C.range, center_coord, stimsize);
    contrast = imresize(squeeze(mean(mean(double(fitmovie(ROIcoord.xvals,ROIcoord.yvals, :))-0.5))), [bins 1],'nearest');
    C_bin = zeros(GLMPars.spikefilters.C.filternumber,bins);
    for i = 1:GLMPars.spikefilters.C.filternumber
       tmp = conv(contrast, C_basis(:,1), 'full');
       C_bin(i,:) = tmp(1:bins); 
    end
end
% NBCoupling 05-28-14
if GLMType.CouplingFilters;
    basis = cp_basis';
    for j_pair=1:GLMPars.spikefilters.cp.n_couplings
        %spikes of neighbor neurons NB
        neighbor_sptimes = neighborspikes.home{j_pair}';
        neighbor_spbins  = ceil(neighbor_sptimes / t_bin);
        neighbor_spbins = neighbor_spbins(find(neighbor_spbins < bins) );
        CP_bin{j_pair}=prep_convolvespikes_basis(neighbor_spbins,basis,bins);
    end
end
if isfield(GLMType, 'Saccades')
   basis = sa_basis';
   spbins = 1:120:bins;
   SA_bin = prep_convolvespikes_basis(spbins, basis, bins);
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end

% Make Covariate Vector of NON stimulus stuff
[paramind] =  prep_paramindGP(GLMType, GLMPars);
p_init     = .01* ones(paramind.paramcount,1);
rawfit.init = p_init;
glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
if isfield(paramind, 'MU')
    glm_covariate_vec( paramind.MU , : ) = MU_bin;
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
if isfield(paramind, 'SA')
    glm_covariate_vec(paramind.SA, :) = SA_bin;
end
if isfield(paramind, 'C')
    glm_covariate_vec(paramind.C, :) = C_bin;
end


% STIM
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,~]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
% clear WN_STA
if GLMType.STA_init && ~strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
    stimsize.width  = size(fitmovie,1);
    stimsize.height = size(fitmovie,2);
    ROIcoord        = ROI_coord(GLMPars.stimfilter.ROI_length, center_coord, stimsize);
    clear stimsize
    STA = WN_STA(ROIcoord.xvals,ROIcoord.yvals,:);
    klen = size(STA,1);
    duration = size(STA, 3);
    STA = reshape(STA, [klen^2,duration])  - mean(STA(:));
    [U,S,V]  = svd (STA);
    imagesc(reshape(U(:,1), [klen, klen]))
    title('Initial Space Filter')
    axis image
    clear X_frame_temp
    p_init(paramind.space1) = U(:,1);
    if strcmp(GLMType.stimfilter_mode, 'rk2')
        p_init(paramind.space2) = U(:,2);
    end
    clear STA U S V
end


% OPTIMIZATION
% INITIALIZE OPTIMIZATION STRUCTURE FOR MATLAB FMIN SOLVERS
if ~GLMType.CONVEX
    GLMPars.optimization.tolfun = GLMPars.optimization.tolfun - 1;
    fittedGLM.GLMPars = GLMPars;
end
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
   'Hessian','off',...
   'largescale','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   );
end
    [pstar fstar eflag output] = fminunc(@(p) glm_SU_optimizationfunction...
        (p,filtertype,paramind,glm_covariate_vec,X_frame,frame_shifts, bpf, home_spbins,t_bin),p_init,optim_struct);

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
end
if isfield(paramind, 'C')
    rawfit.C_basis = C_basis;
    linearfilters.Contrast.Filter     = C_basis * pstar(paramind.C);
    linearfilters.Contrast.startbin   = 1;
    linearfilters.Contrast.note0       = 'Filter starts at "startbin" bins after the spikebin';
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
if isfield(paramind, 'SA')
    rawfit.sa_basis = sa_basis;
    linearfilters.Saccades.Filter = sa_basis * pstar(paramind.SA);
end
% end NBCoupling

% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
center_coord    = glm_cellinfo.slave_centercoord;
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2);
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
if isfield(GLMType, 'postfilter_nonlinearity') && GLMType.postfilter_nonlinearity
    [xvalperformance] = eval_xvalperformance_squared(fittedGLM,testspikes_raster,testmovie,inputstats,neighborspikes.test)
else
    [xvalperformance] = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats,neighborspikes.test)
end
fittedGLM.xvalperformance  = xvalperformance;
 eval(sprintf('save %s/%s.mat fittedGLM',glm_cellinfo.d_save,glm_cellinfo.cell_savename));
printname = sprintf('%s/DiagPlots_%s',glm_cellinfo.d_save,fittedGLM.cellinfo.cell_savename);
printglmfit(fittedGLM,printname)

end

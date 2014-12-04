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
function [fittedGLM] = glm_execute_CP(GLMType, spikes, neighborspikes, fitmovie, glm_cellinfo,troubleshoot)
%% Get rid of all time,  put all inputs into bins

fittedGLM.cellinfo = glm_cellinfo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars           = GLMParams;
fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;
if GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
if GLMType.debug, GLMPars.optimization.tolfun = 3;  end

if GLMType.color
    frames = size(fitmovie, 4);
else
    frames = size(fitmovie,3);
end

bins   = frames * GLMPars.bins_per_frame;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;


% Perhaps we should combine this! With convolving with spikes !
bin_size      = t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end

% NBCoupling
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
clear bin_size basis_params

% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
t_bin        = t_bin;
home_sptimes = spikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );

if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end

% NBCoupling 05-28-14
if GLMType.CouplingFilters;
    n_couplings=length(glm_cellinfo.pairs);
    basis = cp_basis';
    for j_pair=1:n_couplings
        %spikes of neighbor neurons NB
        neighbor_sptimes = neighborspikes.home{j_pair}';
        neighbor_spbins  = ceil(neighbor_sptimes / t_bin);
        neighbor_spbins = neighbor_spbins(find(neighbor_spbins < bins) );
        CP_bin{j_pair}=prep_convolvespikes_basis(neighbor_spbins,basis,bins);
    end
else n_couplings=0;
end
% end NBCoupling

if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end



% Find the correct stimulus related input term
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);

if exist('toubleshoot','var') && troubleshoot.doit
    [X_frame,X_bin]    = prep_stimcelldependentGP(GLMType, GLMPars, fitmovie, center_coord, WN_STA,troubleshoot);
elseif GLMType.color
    [X_frame,X_bin]    = prep_stimcelldependentGP_color(GLMType, GLMPars, fitmovie, center_coord, WN_STA);
else
    [X_frame,X_bin]    = prep_stimcelldependentGP(GLMType, GLMPars, fitmovie, center_coord, WN_STA);
end


clear WN_STA center_coord

%% Set up the parameter index
% NBCoupling
[paramind] =  prep_paramindGP_CP(GLMType, GLMPars, n_couplings); 
p_init     =  .01* ones(paramind.paramcount,1);

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
   'TolX',10^(-(GLMPars.optimization.tolx)));


if GLMType.CONVEX
    glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crashes if not filled out properly
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    X_bin_shift = prep_timeshift(X_bin,shifts);
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
        for j_pair=1:n_couplings
            glm_covariate_vec( paramind.CP{j_pair} , : ) = CP_bin{j_pair};
        end
    end
    % end NBCoupling
    [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),p_init,optim_struct);
end

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
    % NBCoupling 05-28-14
    if isfield(paramind, 'CP')
        for j_pair=1:n_couplings
            convex_cov( paramind.CP{j_pair} , : ) = CP_bin{j_pair};
        end
    end
    filtertype = GLMType.stimfilter_mode;
     %  glm_nonconvex_optimizationfunction(p_init,GLMType.stimfilter_mode,nonconvex_paramind,convex_covariate_vec,X_frame,frame_shifts, bpf, home_spbins,t_bin)
    [pstar fstar eflag output] = fminunc(@(p) glm_nonconvex_optimizationfunction(p,...
        filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins,t_bin),p_init,optim_struct);
end
    


%% Unpack the output into filters
% Do this so we don't have to reinterpret the parameters, just have a final
% filter which we can take home and interpret
% this is admittedly a bit ugly .. 


rawfit.opt_params        = pstar;
rawfit.paramind          = paramind;
rawfit.objective_val     = fstar;

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
        for j_pair=1:n_couplings
            linearfilters.Coupling.Filter{j_pair}     = cp_basis * pstar(paramind.CP{j_pair});
        end
        linearfilters.Coupling.startbin   = 1;  
        linearfilters.Coupling.note = 'Filter starts at "startbin" bins after the spikebin';
    end
    % end NBCoupling


center_coord    = glm_cellinfo.slave_centercoord;
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2); 
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
rawfit.ROIcoord = ROIcoord;
clear stimsize center_coord;
WN_STA           = double(glm_cellinfo.WN_STA); 

if GLMType.color
    for i_filter=1:3
        [STA_sp{i_filter},~]= spatialfilterfromSTA(squeeze(WN_STA(:,:,i_filter,:)),ROIcoord.yvals,ROIcoord.xvals);
    end
else
    [STA_sp,~]= spatialfilterfromSTA(WN_STA,ROIcoord.xvals,ROIcoord.yvals);
end

if GLMType.CONVEX
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')    
        timefilter           = pstar(paramind.X);
        if GLMType.color
            stimfilter = zeros(ROI_length^2,length(paramind.X),3);
            for i_filter=1:3
                stimfilter(:,:,i_filter) = STA_sp{i_filter} * (timefilter');
                linearfilters.Stimulus.space_rk1{i_filter}          = reshape(STA_sp{i_filter}, [ROI_length,ROI_length]);
            end
            stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.X),3]);
            rawfit.spatialfilter = STA_sp;
            linearfilters.Stimulus.Filter             = stimfilter;
            linearfilters.Stimulus.Filter_rank        = 1;
        else
            stimfilter           = STA_sp * (timefilter');
            stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.X)]);
            rawfit.spatialfilter = STA_sp;
            linearfilters.Stimulus.Filter             = stimfilter;
            linearfilters.Stimulus.Filter_rank        = 1;
            linearfilters.Stimulus.space_rk1          = reshape(STA_sp, [ROI_length,ROI_length]);
        end
        linearfilters.Stimulus.time_rk1           = pstar(paramind.X);
        %linearfilters.Stimulus.WN_note            = 'use WN STA as a reference to compare to fitted filters'
        %linearfilters.Stimulus.WN_STA             = WN_STA;
        %linearfilters.Stimulus.WN_STA_space_rk1   = reshape(STA_sp, [ROI_length,ROI_length]);
        %linearfilters.Stimulus.WN_STA_time_rk1    = STA_time;
        linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
        linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
        linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
        linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];
        linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
        linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
        linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
    end
end
if ~GLMType.CONVEX    
    if strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2')       
        timefilter1  = pstar(paramind.time1);
        spacefilter1 = pstar(paramind.space1);
        stimfilter   = spacefilter1 * timefilter1';
        
        if strcmp(GLMType.stimfilter_mode, 'rk2')
            timefilter2  = pstar(paramind.time2);
            spacefilter2 = pstar(paramind.space2);
            stimfilter  = spacefilter1 * timefilter1' + spacefilter2 * timefilter2';
        end
        
        stimfilter = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.time1)]);
        linearfilters.Stimulus.Filter             = stimfilter;
        linearfilters.Stimulus.Filter_rank        = 1;
        linearfilters.Stimulus.time_rk1           = timefilter1;
        linearfilters.Stimulus.space_rk1          = reshape(spacefilter1,[ROI_length,ROI_length]);
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


%% Rudimentary Plots
if exist('troubleshoot','var') && troubleshoot.doit
    clf
    longtitle = sprintf('%s: %s: cid %d',  fittedGLM.cellinfo.exp_nm, fittedGLM.cellinfo.celltype, fittedGLM.cellinfo.cid);

    subplot(3,1,1);
	set(gca, 'fontsize', 10); axis off
	c = 0;
	c=c+1; text(-.1, 1-0.1*c,sprintf('Hacked outputplot: %s',troubleshoot.name));
	c=c+1; text(-.1, 1-0.1*c,longtitle);
	c=c+1; text(-.1, 1-0.1*c,fittedGLM.cellinfo.fitname);
	if isfield(paramind, 'MU')
            c=c+1; text(-.1, 1-0.1*c,sprintf('Tonic Drive: %d', linearfilters.TonicDrive.Filter) );
    end
	c=c+1; text(-.1, 1-0.1*c,sprintf('Fit Date %s',fittedGLM.fit_time));
	c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', fittedGLM.writingcode) );
        
    if isfield(paramind, 'PS')
        subplot(3,2,5)
        plot(linearfilters.PostSpike.Filter)
        title('Post Spike log cif term')
    end
	if isfield(paramind, 'CP')
        subplot(3,2,6)
        % NBCoupling
        plot(linearfilters.Coupling.Filter{j_pair})
        str=[neighborspikes.pair_savename{j_pair}];
        title(str)
        % NBCoupling
    end
        
	stimfilt = linearfilters.Stimulus.Filter;
	subplot(3,2,4);
	imagesc(squeeze(mean(stimfilt,3))); title('stimspacefilter')
	subplot(3,2,3);
	a = squeeze(mean(stimfilt,1));
	timecomponent = mean(a,1);
	plot(timecomponent); title('stimtimefilter')
    orient landscape
    eval(sprintf('print -dpdf %s/%s_glmouputcheck_%s_%s.pdf', troubleshoot.plotdir, troubleshoot.name, fittedGLM.cellinfo.exp_nm, fittedGLM.cellinfo.cell_savename));
    
end
    
    






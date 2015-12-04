% AKHeitman 2015-06-29
% Just a couple changes to ensure negative values.



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

%{
clear; close all
%load('dbug_glmexecute_constrainPS_debugWNfit.mat')
load('dbug_glmexecute_constrainPS_fullWNfit.mat')
domainconstrain_name = 'PS_inhibitorydomainconstrain_post10msec'; 
[fittedGLM] = glm_execute_domainconstrainPS(domainconstrain_name, GLMType,...
    fitspikes_concat,fitmovie_concat,testspikes_raster,testmovie,...
    inputstats,glm_cellinfo); 
%}
function [fittedGLM] = glm_execute_domainconstrainPS(domainconstrain_name, GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,troubleshoot)


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

% IF RUNNING OFF OF INPUT NL LOAD THOSE COEFFS
if isfield(glm_cellinfo, 'GLMPars')
   GLMPars = glm_cellinfo.GLMPars; 
end



fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;
if isfield(GLMType, 'debug') && GLMType.debug
    GLMPars.optimization.tolfun = 1; 
end

frames = size(fitmovie,3);
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
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
clear bin_size basis_params


% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.CouplingFilters;
    basis = cp_basis';
    display('figure out coupling here!  CP_bin');
end

if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


% ORGANIZE STIMULUS COVARIATES
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord

%
if ~GLMType.CONVEX
    GLMPars.optimization.tolfun = GLMPars.optimization.tolfun - 1;
    fittedGLM.GLMPars = GLMPars;
end

%%
% PREPARE PARAMETERS

[paramind] =  prep_paramindGP(GLMType, GLMPars);

% Super Hack to determine PS filter type
% AKHeitman 2015-06-29
if strcmp(domainconstrain_name, 'PS_inhibitorydomainconstrain_post10msec')
    if strcmp(GLMType.fitname_preconstrainPS(end-14:end), '/standardparams')
         lowerbound = -Inf(paramind.paramcount,1);
         upperbound  = Inf(paramind.paramcount,1);
         upperbound(paramind.PS(4:end)) = 0;
    end
    
    fittedGLM.constrained_serach.note = 'how the parameter search was limited in fmincon';
    fittedGLM.constrained_search.lowerbound = lowerbound;
    fittedGLM.constrained_search.upperbound = upperbound;
    
    optim_struct = optimset(...
   'Algorithm','trust-region-reflective',...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   ) ;
end

if strcmp(domainconstrain_name,'PS_netinhibitory_domainconstrain')
    A = zeros(1,paramind.paramcount);
    A(paramind.PS) = sum(ps_basis,1);
    b              = 0;
    
    fittedGLM.constrained_serach.note = 'how the parameter search was limited in fmincon';
    fittedGLM.constrained_search.linearinequality_matrix = A;
    fittedGLM.constrained_search.linearinequality_bound  = b;
    
    
    % Try to figure out how to supply the Hessian
    
    
    %%% OLD VERSION IS REALLY SLOW FOR NSEM %%%'
    % 30 minutes o 3 hours 
    %{
    optim_struct = optimset(...
   'Algorithm','interior-point',...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','fin-diff-grads',...
   'SubproblemAlgorithm', 'cg',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   );
    display('!!!! INTERIOR POINT, FINITE DIFFERENCE OF GRADS!!!')
    %}
end

if strcmp(domainconstrain_name,'PS_netinhibitory_domainconstrain_COB') 
    ps_basis_0 = ps_basis; clear ps_basis
    v        = sum(ps_basis_0,1);
    v        = v / norm(v) ;
    orthog_v = null(v);
    COB      = [v', orthog_v] ;
    ps_basis = (inv(COB) * ps_basis_0')' ;
    
    %%%    
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
    
    
    lowerbound = -Inf(paramind.paramcount,1);
    upperbound  = Inf(paramind.paramcount,1);
    upperbound(paramind.PS(1)) = 0;
    
    fittedGLM.constrained_serach.note = 'how the parameter search was limited in fmincon';
    fittedGLM.constrained_search.lowerbound = lowerbound;
    fittedGLM.constrained_search.upperbound = upperbound;
    
    
    %%%
    optim_struct = optimset(...
   'Algorithm','trust-region-reflective',...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   ) ;
end







%p_init    =  zeros(paramind.paramcount,1);  
p_init     = .01* ones(paramind.paramcount,1);
if isfield(glm_cellinfo, 'p_init')
    p_init = glm_cellinfo.p_init;
    
    if strcmp(domainconstrain_name,'PS_netinhibitory_domainconstrain_COB')
        p_init_psbase = p_init(paramind.PS);
        p_init_new    = inv(COB) * p_init_psbase;        
        p_init(paramind.PS) = p_init_new;
        display('modulating p init to new basis')
    end
    
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
    if isfield(paramind, 'CP')
        glm_covariate_vec( paramind.CP , : ) = CP_bin;
    end
    
    if strcmp(domainconstrain_name, 'PS_inhibitorydomainconstrain_post10msec')...
            || strcmp(domainconstrain_name,'PS_netinhibitory_domainconstrain_COB')
        [pstar fstar eflag output] = fmincon(@(p) ...
            glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),...
            p_init,[],[],[],[],lowerbound,upperbound,[],optim_struct);
    elseif strcmp(domainconstrain_name,'PS_netinhibitory_domainconstrain')
         [pstar fstar eflag output] = fmincon(@(p) ...
            glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),...
            p_init,A,b,[],[],[],[],[],optim_struct);
    end
    [f grad Hess log_cif] = glm_convex_optimizationfunction(pstar,glm_covariate_vec,home_spbins,t_bin);
    
    
    
    %{
    [pstar fstar eflag output] = fminunc(@(p) ...
        glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),...
        p_init,optim_struct);
    [f grad Hess log_cif] = glm_convex_optimizationfunction(pstar,glm_covariate_vec,home_spbins,t_bin);
    plot(ps_basis * pstar_con(paramind.PS),'r'); hold on
    plot(ps_basis * pstar(paramind.PS),'k'); hold on
    % Before constrains
    [pstar fstar eflag output] = fminunc(@(p) ...
        glm_convex_optimizationfunction_constrainPS(p,glm_covariate_vec,home_spbins,t_bin, paramind.PS,ps_basis,PS_balance),...
        p_init,optim_struct);
    [f grad Hess log_cif ps_filter] = glm_convex_optimizationfunction_constrainPS(pstar,glm_covariate_vec,home_spbins,t_bin, paramind.PS,ps_basis,PS_balance);
    %}
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
if isfield(paramind, 'CP')
        rawfit.cp_basis = cp_basis;
        error('need to fill in coupling..  Nice way to handle it')
end

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
[xvalperformance] = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats);
fittedGLM.xvalperformance  = xvalperformance; 
eval(sprintf('save %s/%s.mat fittedGLM',glm_cellinfo.d_save,glm_cellinfo.cell_savename));

currentdir = pwd;
cd(glm_cellinfo.d_save)
printname = sprintf('DiagPlots_%s',fittedGLM.cellinfo.cell_savename);
printglmfit(fittedGLM,printname)
cd(currentdir)

end
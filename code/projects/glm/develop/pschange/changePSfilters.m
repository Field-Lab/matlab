% Call
%{
restoredefaultpath
clear
load figureoutPS.mat

bases_index = [1,2];
changePSfilters(ps_bases,bases_index,GLMType,fitspikes_concat,fitmovie_concat,testspikes_raster,testmovie,inputstats,glm_cellinfo)
%}
function changePSfilters(ps_bases,bases_index,GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo)

%% Setup Covariates
for i_bin = bases_index
    ps_basis = ps_bases{i_bin}.basis;
    ps_name  = ps_bases{i_bin}.note;    

    fittedGLM.cell_savename = glm_cellinfo.cell_savename;
    fittedGLM.d_save        = glm_cellinfo.d_save;
    fittedGLM.cellinfo      = glm_cellinfo;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load up GLMParams compute some universal params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GLMPars           = subR_GLMParams;
    if isfield(GLMType, 'specialchange') && GLMType.specialchange
        GLMPars = subR_GLMParams(GLMType.specialchange_name);
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
    clear bin_size basis_params

    t_bin        = t_bin;
    home_sptimes = fitspikes.home';
    home_spbins  = ceil(home_sptimes / t_bin);
    home_spbins = home_spbins(find(home_spbins < bins) );
    if GLMType.PostSpikeFilter
        basis         = ps_basis';
        PS_bin        = subR_prep_convolvespikes_basis(home_spbins,basis,bins);
    end
    if GLMType.CouplingFilters;
        basis = cp_basis';
        display('figure out coupling here!  CP_bin');
    end

    if GLMType.TonicDrive
        MU_bin = ones(1,bins);
    end
    % PREPARE PARAMETERS
    [paramind] =  subR_prep_paramindGP(GLMType, GLMPars, ps_basis); 
    %p_init     =  zeros(paramind.paramcount,1);  
   


    %% Other Setup stimulus optim struct etc.

    % ORGANIZE STIMULUS COVARIATES
    center_coord       = glm_cellinfo.slave_centercoord;
    WN_STA             = double(glm_cellinfo.WN_STA);
    [X_frame,X_bin]    = subR_prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
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
    p_init     = .01* ones(paramind.paramcount,1);
    if GLMType.CONVEX
        glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
        % Maybe move this inside of the stimulus preparation % 
        bpf         = GLMPars.bins_per_frame;
        shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
        if isfield(GLMPars.stimfilter,'frames_negative')
            shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
        end
        X_bin_shift = subR_prep_timeshift(X_bin,shifts);

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
        if isfield(paramind, 'CP')
            convex_cov( paramind.CP , : ) = CP_bin;
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

    [xvalperformance] = eval_xvalperformance(fittedGLM,testspikes_raster,testmovie,inputstats)
    fittedGLM.xvalperformance  = xvalperformance; 
    eval(sprintf('save %s/%s_%s.mat fittedGLM',glm_cellinfo.d_save,ps_name,glm_cellinfo.cell_savename));
    printname = sprintf('%s/%s_DiagPlots_%s',glm_cellinfo.d_save,ps_name,fittedGLM.cellinfo.cell_savename);
    printglmfit(fittedGLM,printname)
end
end

function GLMPars = subR_GLMParams(optional_change)


%%% First Set all of the default Paramaters
GLMPars.bins_per_frame = 10; 
GLMPars.approx_refresh_hz  = 120;
%GLMPars.dt             = GLMPars.tstim / GLMPars.bins_per_frame;
GLMPars.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMPars.timenotes_1    = 'True tstim is usually .0083275' ;
GLMPars.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMPars.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';


GLMPars.stimfilter.fixedSP_type = 'WNSTA';
GLMPars.stimfilter.ROI_length = 13;  
GLMPars.stimfilter.frames = 30;  
GLMPars.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMPars.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMPars.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMPars.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';


  


GLMPars.optimization.tolfun   = 5;
GLMPars.optimization.tolx     = 9;
GLMPars.optimization.maxiter  = 300; 
GLMPars.optimization.note1    = 'tolfun: significant digits of MLE function being optimized';
GLMPars.optimization.note2    = 'tolx: significant digits of the input variables of the MLE function';

 



%GLMPars.spikefilters.n_psf = 20;  
%GLMPars.spikefilters.n_cp  = 8;
GLMPars.spikefilters.ps_note = 'parameters regarding the post-spike filter';
GLMPars.spikefilters.cp_note = 'parameters regaring coupling filters';
GLMPars.spikefilters.note0 = 'all parameters related to raised sinusoidal humps';
GLMPars.spikefilters.note1 = 'basis built by prep_spikefilterbasisGP / create_histbasis as of 2014-05-3';
GLMPars.spikefilters.ps.ms  = 100 ;      %% post spike filter time length in millisecs
GLMPars.spikefilters.cp.ms  = 100 ;      %% cp spike filter time length in millisecs
%GLMPars.spikefilters.spcng_psf = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
%GLMPars.spikefilters.spcng_cp  = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMPars.spikefilters.BiDirect_CP     = false;
GLMPars.spikefilters.ps.filternumber = 20;
GLMPars.spikefilters.cp.filternumber = 8;
GLMPars.spikefilters.ps.spacing      = pi/2;
GLMPars.spikefilters.cp.spacing      = pi/2;
GLMPars.spikefilters.ps.bstretch     = .05;
GLMPars.spikefilters.ps.alpha        = 0;
GLMPars.spikefilters.cp.bstretch     = .05;
GLMPars.spikefilters.cp.alpha        = 0;
GLMPars.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out



GLMPars.others.fitblockchange = false;


if exist('optional_change','var') && ~isempty(optional_change)
    if strcmp(optional_change, 'ROIlength_9')
        GLMPars.stimfilter.ROI_length = 9;
    elseif strcmp(optional_change, 'BiDirectional_StimFilter')
        GLMPars.stimfilter.frames_negative = GLMPars.stimfilter.frames;
    elseif strcmp(optional_change, 'psfilter_10')
        GLMPars.spikefilters.ps.filternumber = 10;
    elseif strcmp(optional_change, 'Fit_Convergence')
        GLMPars = GLMPars;
    elseif strcmp(optional_change, 'TolFun_7')
        GLMPars.optimization.tolfun = 7;
    elseif strcmp(optional_change, 'TolFun_3')
        GLMPars.optimization.tolfun = 3;
    elseif strcmp(optional_change, 'TolFun_2')
        GLMPars.optimization.tolfun = 2;
    else 
        error('you need to specify how your param changes actually changes the parameters')
    end
end


end
function [spikesconvbasis]  = subR_prep_convolvespikes_basis(binned_spikes,basis,bins)
% AKHeitman 2014-05-04
% Parameter independent!
% basis should be a vector of [basis vectors , bins] already binned
% t_bin is used to put spike times into their appropriate bin 
% t_bin is duration of each bin in msecs
% bins  is the maximum bin number


vectors = size(basis,1); vectorbins= size(basis,2);
%offset by 1,  so PS filter starts in the next bin!!! 
binned_spikes = binned_spikes(find(binned_spikes < bins) ) + 1;


convolvedspikes_base                  = zeros(1,bins);
convolvedspikes_base(binned_spikes+1) = 1; 

convolvedspikes = zeros(vectors, bins + vectorbins - 1);
for i_vec = 1:vectors
    convolvedspikes(i_vec, :) = conv(convolvedspikes_base, basis(i_vec,:), 'full');
end
convolvedspikes = convolvedspikes(:, 1:bins);    

spikesconvbasis = convolvedspikes;
end
function [paramind ] =  subR_prep_paramindGP(GLMType, GLMPars, PS_Basis, CP_Basis)

numParams = 0;


% Assign Parameter Indices.. First to elements that won't brek convexity
% IE. Tonic Drive, Spike filters.. Whatever
if GLMType.TonicDrive
	paramind.MU = 1;
	numParams = numParams+1;
end

if exist('PS_Basis','var')
    PSstart = numParams + 1;
    PSend   = numParams + size(PS_Basis,2);
    paramind.PS = [PSstart:PSend];
    numParams = numParams + length(paramind.PS);
end


% Assign indices to the stim filter
% If CONVEX assign here
if GLMType.CONVEX
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        Xstart     = numParams + 1;  
        Xend       = numParams +  GLMPars.stimfilter.frames;
        paramind.X = [Xstart:Xend];
        numParams  = numParams + GLMPars.stimfilter.frames; 
        if isfield(GLMPars.stimfilter,'frames_negative')
            Xstart2     = numParams + 1;  
            Xend        = numParams +  GLMPars.stimfilter.frames_negative;
            paramind.X  = [Xstart:Xend];
            numParams   = numParams + GLMPars.stimfilter.frames_negative;
        end
    elseif strcmp(GLMType.stimfilter_mode, 'fixedSP-ConductanceBased') 
        paramind.Xnote1 = 'Two Filters but fixed spatial filters, excitatory (STA,time1) and inhibitory (STA,time2)';
        Xstart_1 = numParams + 1;  
        Xend_1   = numParams + GLMPars.stimfilter.frames;
        paramind.time1  = [Xstart_1:Xend_1];
        
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + GLMPars.stimfilter.frames;  
        paramind.time2 = [Xstart_2:Xend_2];
        numParams = numParams  + 2*GLMPars.stimfilter.frames;

        paramind.excitatoryfilter_index = paramind.time1;
        paramind.inhibitoryfilter_index = paramind.time2;
        
        paramind.X = union(paramind.time1,paramind.time2);
    else
        error('you need to properly specifiy the stimfilter in prep_paramind')
    end
end

% If NOT Convex .. Assign here 
if ~GLMType.CONVEX
    convParams = numParams;
    paramind.Xnote0 = 'The stimulus filter X, has non-linear components, breaking convexity';
    
    if strcmp(GLMType.stimfilter_mode, 'rk1') 
        paramind.Xnote1 = 'Stim filter is outer product of space1 and time1';
        Xstart = convParams + 1;  
        Xend   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.X      = [Xstart:Xend];
        paramind.space1 = [Xstart: ((Xstart-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart + GLMPars.stimfilter.ROI_length^2) : Xend ];
        numParams       = convParams  +  GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames;
    end
    if strcmp(GLMType.stimfilter_mode, 'rk2')
        paramind.Xnote1 = 'Stim filter is sum of outer products of (space1,time1) and (space2,time2)';
        Xstart_1 = convParams + 1;  
        Xend_1   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.space1 = [Xstart_1: ((Xstart_1-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart_1 + GLMPars.stimfilter.ROI_length^2) : Xend_1 ];
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;  
        paramind.space2 = [Xstart_2: ((Xstart_2-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time2  = [(Xstart_2 + GLMPars.stimfilter.ROI_length^2) : Xend_2 ];
        numParams = convParams  + 2*( GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames);
    end
    if  strcmp(GLMType.stimfilter_mode, 'rk2-ConductanceBased')
        paramind.Xnote1 = 'Two Filters, excitatory (space1,time1) and inhibitory (space2,time2)';
        Xstart_1 = convParams + 1;  
        Xend_1   = convParams + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;        
        paramind.space1 = [Xstart_1: ((Xstart_1-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time1  = [(Xstart_1 + GLMPars.stimfilter.ROI_length^2) : Xend_1 ];
        Xstart_2 = Xend_1 + 1;  
        Xend_2   = Xend_1 + (GLMPars.stimfilter.ROI_length^2) + GLMPars.stimfilter.frames;  
        paramind.space2 = [Xstart_2: ((Xstart_2-1) + (GLMPars.stimfilter.ROI_length^2))];
        paramind.time2  = [(Xstart_2 + GLMPars.stimfilter.ROI_length^2) : Xend_2 ];
        numParams = convParams  + 2*( GLMPars.stimfilter.ROI_length^2 + GLMPars.stimfilter.frames);

        paramind.excitatoryfilter_index = union(paramind.space1,paramind.time1);
        paramind.inhibitoryfilter_index = union(paramind.space2,paramind.time2);
    end
    
    
    
    paramind.convParams = convParams;
    paramind.convParams_ind = 1:convParams;
end

paramind.paramcount = numParams;

end
function [X_frame,X_bin] = subR_prep_stimcelldependentGPXV(GLMType, GLMPars, stimulus, inputstats, center_coord,STA,troubleshoot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cut the Movie down to ROI
% Normalize Movie and set nullpoint
% output of this section is stim
% stim in [xy, time] coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(stimulus,1);
stimsize.height = size(stimulus,2);
stimsize.frames = size(stimulus,3);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
stim            = stimulus(ROIcoord.xvals, ROIcoord.yvals, :);


fitmoviestats.mean     =  inputstats.mu_avgIperpix;
fitmoviestats.span     =  inputstats.range;
fitmoviestats.normmean =  inputstats.mu_avgIperpix / inputstats.range;

stim   = double(stim);
stim   = stim / fitmoviestats.span;

if strcmp(GLMType.nullpoint, 'mean')
    stim = stim - fitmoviestats.normmean;
else
    error('you need to fill in how to account for stimulus with a different nullpoint')
end




if isfield(GLMType, 'input_pt_nonlinearity') && GLMType.input_pt_nonlinearity
   % display('implementing nonlinearity')
    newstim = stim;
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        par       = GLMPars.others.point_nonlinearity.increment_to_decrement;
        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ;
        
        pos_stim          = find(stim > 0 );
        neg_stim          = find(stim < 0 );
        
        newstim(pos_stim) = pos_mult * (newstim(pos_stim)); 
        newstim(neg_stim) = neg_mult * (newstim(neg_stim));
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
        par = GLMPars.others.point_nonlinearity.increment_to_decrement;
        par_shift = GLMPars.others.point_nonlinearity.shiftmean; 

        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ;
        
        pos_stim          = find(stim > par_shift);
        neg_stim          = find(stim < par_shift);
        
        newstim(pos_stim) = pos_mult * (newstim(pos_stim)-par_shift) + par_shift; 
        newstim(neg_stim) = neg_mult * (newstim(neg_stim)-par_shift) + par_shift;        
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'log')
        display('implenting log')
        newstim = stim + fitmoviestats.normmean + 1; % now back on 0 1 scale
        newstim = log(newstim);
        newstim = newstim - log(fitmoviestats.normmean+1);
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'exp')
        display('implenting exp')
        newstim = stim + fitmoviestats.normmean; % now back on 0 1 scale
        newstim = exp(newstim);
        newstim = newstim - exp(fitmoviestats.normmean); 
    elseif strcmp(GLMType.input_pt_nonlinearity_type,'piecelinear_fourpiece_eightlevels')
        newstim = stim + fitmoviestats.normmean;
        
        quartile_1 = find(newstim<=.25);
        quartile_2 = setdiff(find(newstim<=.5),quartile_1);
        
        quartile_4 = find(newstim>.75);
        quartile_3 = setdiff(find(newstim >.5),quartile_4);
        
        
        coeff = GLMPars.others.point_nonlinearity.coefficients;
        slope1 = coeff.slope_quartile_1;
        slope2 = coeff.slope_quartile_2;
        slope3 = coeff.slope_quartile_3;
        slope4 = coeff.slope_quartile_4;
        
        offset1 = 0;
        offset2 = .25* slope1;
        offset3 = .25*(slope1+slope2);
        offset4 = .25*(slope1+slope2+slope3);
        
        newstim(quartile_1) = slope1 * (newstim(quartile_1) - .00) + offset1;
        newstim(quartile_2) = slope2 * (newstim(quartile_2) - .25) + offset2;
        newstim(quartile_3) = slope3 * (newstim(quartile_3) - .50) + offset3;
        newstim(quartile_4) = slope4 * (newstim(quartile_4) - .75) + offset4;
        
        newstim = newstim - offset2;
    else
        display('error, need to properly specifiy input non-linearity')
    end
    stim = newstim; clear newstim
end



%}
stim = reshape(stim, [ROI_length^2 , stimsize.frames]);

clear stimsize ROIcoord ROI_length fitmoviestats 

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create spfilter from the STA 
% output: spfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear') 
    ROI_length      = GLMPars.stimfilter.ROI_length;
    stimsize.width  = size(stimulus,1);
    stimsize.height = size(stimulus,2); 
    ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
    spfilter = spatialfilterfromSTA(STA,ROIcoord.xvals,ROIcoord.yvals);
    if exist('troubleshoot','var') && troubleshoot.doit
        clf
        subplot(3,2,[1 2]);   set(gca, 'fontsize', 10); axis off
        c = 0;
        c=c+1; text(-.1, 1-0.1*c,sprintf('Trouble shooting: %s',troubleshoot.name));
        c=c+1; text(-.1, 1-0.1*c,sprintf('Specifically: spfilter from WN-STA "glm-nospace/spatialfilterfromSTA"' ));
        c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
        c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')) );


        subplot(3,2,[3 5]);  set(gca, 'fontsize', 10); imagesc(reshape(spfilter,[ROI_length, ROI_length])); colorbar
        xlabel('pixels'); ylabel('pixels');
        title('Spatial Filter first rank of the STA');

        subplot(3,2,[4 6]); set(gca, 'fontsize', 10);   imagesc( (squeeze(mean(STA,1))' )); colorbar
        ylabel('frames');  xlabel('pixel rows')
        title('Raw STA, columns collapsed to 1-spatial dimension');


        orient landscape
        eval(sprintf('print -dpdf %s/%s_prepstimcellGP_spfilterfromSTA.pdf', troubleshoot.plotdir, troubleshoot.name));
    end 
    spfilter  = spfilter';  % dimension [1,ROI_length]
end

clear stimsize ROIcoord ROI_length fitmoviestats 

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Creat the final stim dependent input (X) to the optimization algorithm
% X_bin  ([rank, bins])
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
    if strcmp(GLMPars.stimfilter.fixedSP_type, 'WNSTA')
        X_frame = (spfilter) * stim;
    end
elseif strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2') || ...
        strcmp(GLMType.stimfilter_mode, 'rk2-ConductanceBased')||strcmp(GLMType.stimfilter_mode, 'rk1-newfit')||strcmp(GLMType.stimfilter_mode, 'fullrank')
    X_frame = stim;
else
    error('you need to tell prep_stimcelldependentGP how to process stim for your spatial filter')
end

frames = size(X_frame,2);
dim    = size(X_frame,1);

bpf    = GLMPars.bins_per_frame;
bins   = bpf * frames;
X_bin  = repmat(X_frame, [ bpf,1]); 
X_bin  = reshape(X_bin, [dim, bins]);

clear frames bpf bins 


end
%{
function basis_vectors = subR_prep_spikefilterbasisGP(basis_params,bin_size)
% MY_FUNCTION     This finds all the params associated with the 2 filters.
%                 Whcih were not specified in   
%   Calls create_histbasis  (should be renamed create_psbasis)
%
%                                           'bore' - activate remotely
%  NOT A VERY GNERAL FUNCTION
%
% Works as of 
%



% Done in an microbin that's been adjusted for with the ps_fratio factor
%ps_f ration mean spost spike filter to linear filter ratio!!!
%init_pars.ps_timebins = floor((init_pars.spikebins_perstimframe*init_pars.k_stimframes)*init_pars.ps_fratio); % on fine timescale ~ 6/factor m(s)^-1

dt = bin_size;
 
basis_params.timebins  = floor( (basis_params.ms/1000) / dt );
basis_params.beta    = (basis_params.timebins-1)*dt; % ending point for ps filters
% History filters (postspike and coupling)
bstretch = basis_params.bstretch;
alpha_ps = basis_params.alpha;
beta_ps  = basis_params.beta;
if (basis_params.filternumber > 0)
   basis_vectors = create_histbasis(alpha_ps,beta_ps,bstretch,...
      basis_params.filternumber,basis_params.timebins,dt,'L2',basis_params.spacing);
else
   basis_params.basis = [];
end

end
%}
function X_bin_shift = subR_prep_timeshift(X_bin,shifts)

bins = size(X_bin,2);
dim  = size(X_bin,1);

X_bin_shift = zeros(dim*length(shifts) , bins);

for i_shift = 1:length(shifts)
    step = shifts(i_shift);
    index = (i_shift-1)*dim +  [1:dim] ;
    X_bin_shift(index,:) = circshift(X_bin , [0 step]);
end

end 
function [ROIcoord] = ROI_coord(ROI_length, center, stimsize)
% ROI_length preferably odd
% Stimsize.width (corresponds to x)
% Stimsizse.height  (corresponds to y)
% center.x_coord, center.y_coord


% ROI size cannot exceed screen dimensions
if ROI_length > stimsize.width || ROI_length > stimsize.height
    ROI_length = min(stimsize.width, stimsize.height);
end

% Initial coord
modvec = -floor(ROI_length/2) : floor(ROI_length/2) ;
xdim = center.x_coord + modvec;
ydim = center.y_coord + modvec;

% Correct if ROI falls off the screen
if min(xdim)<1
	xdim = 1:ROI_length;
end
if min(ydim)<1
	ydim = 1:ROI_length;
end
if max(xdim) > stimsize.width
	xdim = (stimsize.width - ROI_length + 1):stimsize.width ;
end
if max(ydim) > stimsize.height
	ydim = (stimsize.height - ROI_length + 1):stimsize.height ;
end

% Assign the final answer
ROIcoord.xvals =  xdim;
ROIcoord.yvals  = ydim;


end
%{
function [basis,t_psb] = create_histbasis(alpha,beta,bstretch,nofilters,Mhist,dt,norm_mode,spacing)
% Directly From CHAITU
%basis =  construct_cosine_basis2(alpha,beta,bstretch,nofilters,dt,spacing); % raised cosine bumps
[basis,~,t_psb] =  subR_construct_cosine_basis3(alpha,beta,bstretch,nofilters,dt,spacing); % edoi, 2011-01-07.

if (alpha < 0)
    offset = floor(abs(alpha)/dt);
    basis1 = basis(offset+1:min(size(basis,1),offset+Mhist),:);
    basis = [basis1; zeros(max(0,Mhist-size(basis1,1)),size(basis,2))];
else
    basis = [zeros((Mhist-size(basis,1)),size(basis,2)); basis];
end


switch(norm_mode)
    case 'area'
        areas = sum(basis,1); % Normalization by signed area
    case 'L2'
        areas = sqrt(sum(basis.^2,1)); % Normalization by L2 norm
    case 'none'
        areas = ones(1,size(basis,2));
end
basis = basis ./ repmat(areas,Mhist,1); % normalize the basis

basis(isnan(basis)) = 0;

end
function [basis, phi,tvec] = subR_construct_cosine_basis3(alpha,beta,b,nfilters,dt,spacing,fl_out)
% Directly From CHAITU
% functional form of raised cosine (with logarithmic time scale) is
% g_j(t) = (1/2)cos(a*log[b*t + c] - phi_j) + (1/2)
% for t such that (phi_j - pi) <= a*log(b*t + c) <= (phi_j + pi)

% alpha - starting point (can be negative)
% beta - finishing point
% b - time scaling INSIDE log
% nfilters - number of filters to construct
% dt - timestep (frame or spike -- note this parameter is working weird..)

% based on construct_cosine_basis2.m, edoi, 2011-01-07.
% add: in the first basis, basis element is set 1 until the peak (as in Pillow et al, 2008, suppl.)


if ~exist('fl_out','var')
   fl_out = 0;
end

tvec = alpha:dt:beta;
% define c s/t [b*t + c] = dt when t = alpha; this is necessary for
% the logarithmic time scaling (where t need to be > 0)
c = -b*alpha + dt;
if (~exist('spacing','var'))
   spacing = pi/2; % default spacing (square sum is constant)
end
a = (spacing*(nfilters-1) + 2*pi) * 1/log((b*beta + c)/(b*alpha + c));
minphi = a*log(b.*alpha + c) + pi;
maxphi = a*log(b.*beta  + c) - pi;

if (maxphi < minphi)
   fprintf('ERROR: alpha beta are unsuitable for cosine basis construction, minphi=%f, maxphi=%f !\n',...
      minphi,maxphi);
   return;
end

T = length(tvec);

phi = linspace(minphi,maxphi,nfilters);
if (length(phi) > 1)
   A = norm(diff(phi) - repmat(spacing,1,nfilters-1));
   if (A > 10^(-6))
      fprintf('ERROR: spacing is incorrect to precision %f!n',A);
   end
end

basis = zeros(T,nfilters);
% construct the postspike basis functions
arg = a.*log(b.*tvec + c);
for k=1:nfilters
   idx_k =  find((arg > phi(k)-pi) & (arg < phi(k)+pi));
   if isempty(idx_k)
      fprintf('timespan of filter %d is empty!\n',k);
   else
      if fl_out == 1
         fprintf('timespan of filter %d is (%f,%f)\n',k,tvec(idx_k(1)),tvec(idx_k(end)));
      end
      basis(idx_k,k) = 0.5.*(cos(arg(idx_k) - phi(k))+1);
      if k == 1 % add
         idx_k0 =  (arg < phi(1));
         basis(idx_k0,1) = 1;
      end
   end
end

% correction of the initial sharp raise: edoi, 2012-01-20
nbin_ar = find(tvec<3/4*10^-3,1,'last'); % for modeling the absolute refractory period
if sum(idx_k0) > nbin_ar
   nbin_shift = sum(idx_k0) - nbin_ar;
   basis = [basis(nbin_shift+1:end,:);zeros(nbin_shift,nfilters)];
end

end
%}
function [leadspfilter leadtimefilter] = spatialfilterfromSTA(STA,xcoord,ycoord)


% Reshape into space,time  2d notation
duration = size(STA,3);
klen = length(xcoord);
STA = (STA(xcoord,ycoord,:) );
STA = reshape(STA, [klen^2,duration])  - mean(STA(:)) ;

% Making sure no wiered NAN stuff
isfiniteSTA = isfinite(STA);

% Singular Value Decomposition
if isempty(find(isfiniteSTA == 0 ) )
    [U,S,V]  = svd (STA);
    S = diag(S);
    % Choosing the V(5,1) lets us make the direction consistent
    % OFFs will look OFF , ONs will look ON
    xx = ( S(1)*U(:,1)*V(5,1) ) / norm( S(1)*U(:,1)*V(5,1) ) ;
    leadspfilter = xx;
    
    yy = S(1)*V(:,1) / norm( S(1)*V(:,1) ); 
    leadtimefilter = yy;
else
    error('STA is not well definted')
end



end
function  [f grad Hess log_cif]= glm_convex_optimizationfunction(linear_params,covariates,spikebins,bin_duration)
 %%% PURPOSE %%%
% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized, column vector
% Covariates: time dependent input with multiplies into the params
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 


% AKHEITMAN 2014-12-04
%%
% Initialize
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;


% Find Conditional Intensity and its log
lcif = p' * COV;
cif  = exp(lcif);


% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);

% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );

% Evaluate the hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');


% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end
function [xvalperformance] = eval_xvalperformance(fittedGLM,testspikes,testmovie,inputstats)
%%
bpf               = fittedGLM.bins_per_frame;
params.bindur     = fittedGLM.t_bin;
params.bins       = fittedGLM.bins_per_frame *size(testmovie,3); 
params.trials     = length(testspikes.home); 
params.frames     = size(testmovie,3);
params.testdur_seconds = params.bindur * params.bins ;   
center_coord = fittedGLM.cellinfo.slave_centercoord;
teststim     = testmovie;
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

%%
logicalspike = zeros(params.trials,params.bins) ;         
for i_blk = 1 : params.trials
    spt = testspikes.home{i_blk};
    binnumber = ceil(spt / params.bindur );
    logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
end 
clear i_blk spt sptimes

%%
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
[X_frame] = subR_prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, teststim,inputstats,center_coord) ;
clear GLMType_fortest
GLMType = fittedGLM.GLMType;


  
    %% Set up CIF Components
    

MU = fittedGLM.linearfilters.TonicDrive.Filter;
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end
K  = fittedGLM.linearfilters.Stimulus.Filter;

% HUGE HACK AKHeitman 2014-10-21
% rk1 filters are misscaled... too hard to dig out
% rk1 filters are fit fine
% this is confirmed to be the correct factor though!
%if strcmp(fittedGLM.GLMType.stimfilter_mode, 'rk1')
%    K = 2*K;
%end
K  = reshape(K, [ROI_pixels, length(frame_shifts)]);



KX = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = subR_prep_timeshift(X_frame(i_pixel,:),frame_shifts);
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
if GLMType.PostSpikeFilter
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
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


%%
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


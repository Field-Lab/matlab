% AKHeitman 2014-04-14
% Called by glm_nospace and glm_space
% Conscise location of detailed parameters of the GLM fit

% bins per frame is literally the number of bins (used to hold spikes) per
% frame

%%%%%% DEFAULT VALUES AS DETERMINED 2013-01-07 !!!!
function GLMPars = GLMParams(optional_change)


%%% First Set all of the default Paramaters
GLMPars.bins_per_frame = 10; %orig 10
GLMPars.approx_refresh_hz  = 120;
%GLMPars.dt             = GLMPars.tstim / GLMPars.bins_per_frame;
GLMPars.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMPars.timenotes_1    = 'True tstim is usually .0083275' ;
GLMPars.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMPars.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';


GLMPars.stimfilter.fixedSP_type = 'WNSTA';

GLMPars.stimfilter.ROI_length = 13;  



GLMPars.stimfilter.frames = 30;  % orig 30
GLMPars.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMPars.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMPars.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMPars.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';
% GLMPars.stimfilter.frames_negative = 2;

% NBsubunits eventually these will be parameters in a dual fitting type thing
GLMPars.others.point_nonlinearity.increment_to_decrement=3;
GLMPars.others.point_nonlinearity.scalar_raisedpower=2;
GLMPars.subunit.size = 3;
GLMPars.subunit.pretime_filter_frames = 30;
%GLMPars.time = 'pre_conv';
% GLMPars.subunit.time_after = 'fit';
  


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
GLMPars.spikefilters.ps.filternumber = 20; %orig 20 %10 with bstretch 0.95 works pretty well
% GLMPars.spikefilters.ps.filternumber = 1;
% GLMPars.spikefilters.cp.filternumber = 4;
GLMPars.spikefilters.cp.filternumber = 8;
GLMPars.spikefilters.ps.spacing      = pi/2;
GLMPars.spikefilters.cp.spacing      = pi/2;
GLMPars.spikefilters.ps.bstretch     = 0.05;
GLMPars.spikefilters.ps.alpha        = 0;
GLMPars.spikefilters.cp.bstretch     = .05;
GLMPars.spikefilters.cp.alpha        = 0;
GLMPars.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.n_couplings = 6  ;

GLMPars.saccadefilter.ms = 100;
GLMPars.saccadefilter.filternumber = 4;
GLMPars.saccadefilter.spacing = pi;
GLMPars.saccadefilter.bstretch = 0.5;
GLMPars.saccadefilter.alpha = 0;
GLMPars.saccadefilter.fratio = 0.95;

GLMPars.spikefilters.C.ms = 500;
GLMPars.spikefilters.C.filternumber = 5;
GLMPars.spikefilters.C.spacing = pi/2;
GLMPars.spikefilters.C.bstretch = 0.75;
GLMPars.spikefilters.C.alpha = 0;
GLMPars.spikefilters.C.fratio = 0.5;
GLMPars.spikefilters.C.range = 20;

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

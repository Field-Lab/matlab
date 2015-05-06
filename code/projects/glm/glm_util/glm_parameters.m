
%% GLM Type

% Spatial Filter Type. Choose FixedSP, Rk1 or Rk2
% GLMType.stimfilter_mode = 'fixedSP_rk1_linear'; GLMType.CONVEX = true;
% GLMType.stimfilter_mode = 'rk1'; GLMType.CONVEX = false; 
GLMType.stimfilter_mode = 'rk2'; GLMType.CONVEX = false;

% Coupling on or off?
GLMType.CouplingFilters = false;

% You probably don't want to change these.
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';
GLMType.TonicDrive = true; % background firing
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.nullpoint  = 'mean';
GLMType.map_type   = 'mapPRJ';
GLMType.debug      = false;
GLMType.specialchange = false;

%% GLM PARAMETERS
% Detailed parameters of the model, including time bins, filter bases, etc.
% You probably don't want to mess with these much, but if you want to tweak
% the model fits, here ya go.

% Timing
GLMPars.bins_per_frame = 10; %orig 10
GLMPars.approx_refresh_hz  = 120;
%GLMPars.dt             = GLMPars.tstim / GLMPars.bins_per_frame;
GLMPars.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMPars.timenotes_1    = 'True tstim is usually .0083275' ;
GLMPars.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMPars.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';

% STA and Spatial Filter Size
GLMPars.stimfilter.fixedSP_type = 'WNSTA';
GLMPars.stimfilter.ROI_length = 5;  
GLMPars.stimfilter.frames = 30;  % orig 30
GLMPars.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMPars.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMPars.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMPars.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';

% Optimization
GLMPars.optimization.tolfun   = 5;
GLMPars.optimization.tolx     = 9;
GLMPars.optimization.maxiter  = 300; 
GLMPars.optimization.note1    = 'tolfun: significant digits of MLE function being optimized';
GLMPars.optimization.note2    = 'tolx: significant digits of the input variables of the MLE function';

% Post-spike and Coupling Basis settings
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
GLMPars.spikefilters.ps.bstretch     = 0.05;
GLMPars.spikefilters.ps.alpha        = 0;
GLMPars.spikefilters.cp.bstretch     = 0.05;
GLMPars.spikefilters.cp.alpha        = 0;
GLMPars.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.n_couplings = 6  ;

GLMPars.others.fitblockchange = false;


function [GLMT, GLMP] = glm_parameters
%% GLM Type

% Spatial Filter Type. Choose FixedSP, Rk1 or Rk2

% GLMT.stimfilter_mode = 'fixedSP_rk1_linear'; GLMT.CONVEX = true;
GLMT.stimfilter_mode = 'rk1'; GLMT.CONVEX = false; 
% GLMT.stimfilter_mode = 'rk2'; GLMT.CONVEX = false;

% Coupling on or off?
GLMT.CouplingFilters = true;

% You probably don't want to change these.
GLMT.cone_model = '8pix_Identity_8pix'; GLMT.cone_sname='p8IDp8';
GLMT.TonicDrive = true; % background firing
GLMT.StimFilter = true;
GLMT.PostSpikeFilter = true;
GLMT.nullpoint  = 'mean';
GLMT.map_type   = 'mapPRJ';
GLMT.debug      = false;
GLMT.specialchange = false;

%% GLM PARAMETERS
% Detailed parameters of the model, including time bins, filter bases, etc.
% You probably don't want to mess with these much, but if you want to tweak
% the model fits, here ya go.

% Timing
GLMP.bins_per_frame = 10; %orig 10
GLMP.approx_refresh_hz  = 120;
%GLMP.dt             = GLMPars.tstim / GLMPars.bins_per_frame;
GLMP.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMP.timenotes_1    = 'True tstim is usually .0083275' ;
GLMP.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMP.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';

% STA and Spatial Filter Size
GLMP.stimfilter.fixedSP_type = 'WNSTA';
GLMP.stimfilter.ROI_length = 13;  
GLMP.stimfilter.frames = 30;  % orig 30
GLMP.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMP.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMP.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMP.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';

% Optimization
GLMP.optimization.tolfun   = 5;
GLMP.optimization.tolx     = 9;
GLMP.optimization.maxiter  = 300; 
GLMP.optimization.note1    = 'tolfun: significant digits of MLE function being optimized';
GLMP.optimization.note2    = 'tolx: significant digits of the input variables of the MLE function';

% Post-spike and Coupling Basis settings
GLMP.spikefilters.ps_note = 'parameters regarding the post-spike filter';
GLMP.spikefilters.cp_note = 'parameters regaring coupling filters';
GLMP.spikefilters.note0 = 'all parameters related to raised sinusoidal humps';
GLMP.spikefilters.note1 = 'basis built by prep_spikefilterbasisGP / create_histbasis as of 2014-05-3';
GLMP.spikefilters.ps.ms  = 100 ;      %% post spike filter time length in millisecs
GLMP.spikefilters.cp.ms  = 100 ;      %% cp spike filter time length in millisecs

%GLMP.spikefilters.spcng_psf = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
%GLMP.spikefilters.spcng_cp  = pi/2;  %% it could be set as pi, but pi/2 is better for "uniform" sampling.
GLMP.spikefilters.BiDirect_CP     = false;
GLMP.spikefilters.ps.filternumber = 10;
GLMP.spikefilters.cp.filternumber = 8;
GLMP.spikefilters.ps.spacing      = pi/2;
GLMP.spikefilters.cp.spacing      = pi/2;
GLMP.spikefilters.ps.bstretch     = 0.05;
GLMP.spikefilters.ps.alpha        = 0;
GLMP.spikefilters.cp.bstretch     = 0.05;
GLMP.spikefilters.cp.alpha        = 0;
GLMP.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMP.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out
GLMP.spikefilters.cp.n_couplings = 6  ;

GLMP.others.fitblockchange = false;

end


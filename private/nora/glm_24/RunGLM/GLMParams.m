% AKHeitman 2014-04-14
% Called by glm_nospace and glm_space
% Conscise location of detailed parameters of the GLM fit

% bins per frame is literally the number of bins (used to hold spikes) per
% frame

%%%%%% DEFAULT VALUES AS DETERMINED 2013-01-07 !!!!
function GLMPars = GLMParams(optional_change)


%%% First Set all of the default Paramaters
GLMPars.bins_per_frame = 10; 
GLMPars.approx_refresh_hz  = 120;
%GLMPars.dt             = GLMPars.tstim / GLMPars.bins_per_frame;
GLMPars.timenotes_0    = 'tstim is ~time in seconds of frame refresh,  dt is the ~time per GLM bin';
GLMPars.timenotes_1    = 'True tstim is usually .0083275' ;
GLMPars.timenotes_2    = 'true tstim is measured from the triggers in each datarun, usually by the Directories_Params function';
GLMPars.timenotes_3    = 'true tstim only matters for binning the spike times when we exceed timescales of seconds';


GLMPars.stimfilter.fixedSP_type = 'WNSTA';
GLMPars.stimfilter.ROI_length = 11; %originally 13
GLMPars.stimfilter.frames = 30;
GLMPars.stimfilter.note1 = 'ROI_length: refers to dimension of stimulus used for GLM fitting';
GLMPars.stimfilter.note2 = 'ROI_length: will also be size of spatial filter if we are fitting a spatial filter';
GLMPars.stimfilter.note3 = 'Frames: Time duration of the fitted stim filter in frames';
GLMPars.stimfilter.note4 = 'Frames: Time duration of the fitted stim filter in frames';

% NBsubunits eventually these will be parameters in a dual fitting type thing
GLMPars.others.point_nonlinearity.increment_to_decrement=3;
GLMPars.others.point_nonlinearity.scalar_raisedpower=2;
GLMPars.subunits=ones(2,2);

% optimization
GLMPars.optimization.tolfun = 5;
GLMPars.optimization.tolx   = 9;
GLMPars.optimization.maxiter = 100; 
GLMPars.optimization.note1 = 'tolfun: significant digits of MLE function being optimized';
GLMPars.optimization.note2 = 'tolx: significant digits of the input variables of the MLE function';


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
GLMPars.spikefilters.BiDirect_CP = false;
GLMPars.spikefilters.ps.filternumber = 20; %20
GLMPars.spikefilters.cp.filternumber = 4; %8
GLMPars.spikefilters.ps.spacing      = pi/2; %pi/2
GLMPars.spikefilters.cp.spacing      = pi/2; %pi/2
GLMPars.spikefilters.ps.bstretch     = 0.05; %0.05
GLMPars.spikefilters.ps.alpha        = 0; %0
GLMPars.spikefilters.cp.bstretch     = 0.05; %0.05
GLMPars.spikefilters.cp.alpha        = 0; %0
GLMPars.spikefilters.ps.fratio = .5  ;  % legacy afraid to take out
GLMPars.spikefilters.cp.fratio = .4  ;  % legacy afraid to take out

GLMPars.others.fitblockchange = false;

if exist('optional_change','var') && ~isempty(optional_change)
    if strcmp(optional_change, 'ROIlength_9')
        GLMPars.stimfilter.ROI_length = 9;
    elseif strcmp(optional_change, 'extra_coupling')
    else 
        error('you need to specify how your param changes actually changes the parameters')
    end
end


end

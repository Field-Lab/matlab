function datarun = cone_preprocess(datarun, varargin)
% CONE_PREPROCESS   Put stimulus into cone space and bin spikes into stimulus frame epochs
% usage: datarun = cone_preprocess(datarun, opts)
%
% DATARUN must include specific information in DATARUN.STIMULUS:
%   DATARUN.STIMULUS.REFRESH_TIME       The refresh period of the stimulus in ms
%   DATARUN.STIMULUS.STIMOBJECT         The object to get frames from, e.g. JAVA_MOVIE or a raw array
%   DATARUN.STIMULUS.GET_MOVIE_FRAME    Function handle for getting frame I from STIMOBJECT; see 
%                                       CONEPREPROCESS_WNM.GET_MOVIE_FRAME_FROM_JAVA_MOVIE for example.
%
% Additionally, we may require:
%   DATARUN.STIMULUS.DURATION       Used to set end_stim if no OPTS.END_TIME given.
%   DATARUN.STIMULUS.WC_PREPROCESS  Function handle for preprocessing the cone WC data; it may be more 
%                                   efficient to preprocess this than to do additional processing on 
%                                   every movie frame; see CONEPREPROCESS_WNM.WC_PREPROCESS_FOR_JAVA_MOVIE
%                                   for example.  Basically, needs to complement GET_MOVIE_FRAME logic.
%   DATARUN.STIMULUS.FIELD_WIDTH    May be needed for WC_PREPROCESS
%   DATARUN.STIMULUS.FIELD_HEIGHT   May be needed for WC_PREPROCESS
%
% For repeat stimuli, provide:
%   DATARUN.STIMULUS.TRIGGERS_PER_TRIAL    Hopefully there are exactly the same number of triggers per
%                                          repeat trial, otherwise this will get hairy...
%
% Additional OPTS:
%   cone_data_ind   []                  Passed to LOAD_CONES if necessary; needed if there are multiple 
%                                       cone findings saved for this run
%   start_time      0                   Start time to process, in seconds
%   end_time        DATARUN.DURATION    End time, in seconds
%
% See also: CALC_CONE_INPUTS, GET_BINNED_SPIKES, CONEPREPROCESS_SAVE, CONEPREPROCESS_WNM
%
% 2012-09, phli
%

opts = inputParser();
opts.addParamValue('conepath', []);
opts.addParamValue('cone_data_ind', []);
opts.addParamValue('start_time', 0);
opts.addParamValue('end_time', datarun.stimulus.duration);
opts.parse(varargin{:});
opts = opts.Results;


% Load data as necessary
if ~isfield(datarun, 'triggers') || ~isfield(datarun, 'spikes') || ~isfield(datarun, 'duration')
    datarun = load_neurons(datarun);
end

% Calculate the desired stim frames
start_stim = floor(1+opts.start_time/datarun.stimulus.refresh_time);
end_stim = floor(1+opts.end_time/datarun.stimulus.refresh_time);
stims = start_stim:end_stim;

% load cone info
if isempty(opts.conepath)
    if ~isfield(datarun, 'cones') || isempty(datarun.cones.centers)
        datarun = load_cones(datarun, opts.cone_data_ind);
    end
    
    path2data=datarun.names.rrs_prefix;
    tmp=regexp(path2data,'/');
    tmp=path2data(1:tmp(4));    
    
%     tmp=regexp(datarun.names.rrs_prefix,'data');
%     tmp=datarun.names.rrs_prefix(1:tmp-1);
    
    opts.conepath = datarun.names.cones;
    opts.conepath = fullfile(tmp, opts.conepath);
end
% opts.conepath = fullfile(single_cone_path(), opts.conepath);

if ~exist('Wc','var')
    load(fullfile(opts.conepath, 'Wc.mat'));
end
 
% Preprocess Wc for efficiency if desired
if isfield(datarun.stimulus, 'wc_preprocess') && ~isempty(datarun.stimulus.wc_preprocess)
    Wc = datarun.stimulus.wc_preprocess(Wc, datarun.stimulus);
end

% Reexpress movie in terms of cone signals
datarun.cone_inputs = calc_cone_inputs(Wc, stims, datarun.stimulus.stimobject, datarun.stimulus.get_movie_frame);

% get the spikes in each bin
if isfield(datarun.stimulus, 'triggers_per_trial') && ~isempty(datarun.stimulus.triggers_per_trial)

		
    triggers      = datarun.triggers;
    trigspertrial = datarun.stimulus.triggers_per_trial;

    % Check whether trigspertrial divides evenly into triggers
    rounded = floor(length(triggers) / trigspertrial) * trigspertrial;
    if rounded < length(triggers)
        warning('Triggers per trial did not divide evenly into triggers; was the run interrupted, or did you use the wrong triggers per trial?');
        triggers = triggers(1:rounded);
    end

    triggers = reshape(triggers, trigspertrial, []);
    ntrials = size(triggers,2);

    spikes = cell(ntrials,1);
    for i = 1:ntrials
        temp = get_binned_spikes(datarun, stims, 'triggers', triggers(:,i));
        spikes{i} = temp.spike_rate;
    end
    
    datarun.spike_rate = spikes;

else
    datarun = get_binned_spikes(datarun, stims);
end

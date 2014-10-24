function datarun = conepreprocess_wnrepeats(datarun, moviexmlpath, triggers_per_trial, varargin)
% CONEPREPROCESS_WNREPEATS    Setup datarun.stimulus for CONE_PREPROCESS on white noise repeats
% usage: datarun = conepreprocess_wnrepeats(datarun, moviexmlpath, triggers_per_trial, opts)
%
% Sets up DATARUN.STIMULUS by loading the java_movie and then setting a 
% few additional fields.  OPTS are simply passed on to CONE_PREPROCESS.
%
% See also: CONEPREPROCESS_WNM, CONE_PREPROCESS, CALC_CONE_INPUTS, GET_BINNED_SPIKES
%
% 2012-09 phli
%

if ~isfield(datarun, 'stimulus') || ~isfield(datarun.stimulus, 'java_movie')
    datarun = load_java_movie(datarun, moviexmlpath, datarun.triggers(1:triggers_per_trial));
    
    % FIXME: This version of load_java_movie does shitty hardcoded things
    datarun.stimulus = rmfield(datarun.stimulus, {'interval', 'monitor_refresh'});
end
datarun.stimulus.refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
datarun.stimulus.stimobject = datarun.stimulus.java_movie;
datarun.stimulus.get_movie_frame = @get_movie_frame_from_java_movie;
datarun.stimulus.wc_preprocess = @wc_preprocess_javamovie_rawmovie;

% FIXME: Not quite as long as it could be; need to get the frames setting
% from stim lisp to know how long the stim was really still running
datarun.stimulus.duration = diff(datarun.triggers([1 triggers_per_trial]));
datarun.stimulus.triggers_per_trial = triggers_per_trial;

datarun = cone_preprocess(datarun, varargin{:});


function frame = get_movie_frame_from_java_movie(java_movie, i)
frame = double(java_movie.getFrame(i-1).getBuffer());
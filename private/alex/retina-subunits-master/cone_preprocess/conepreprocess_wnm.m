function datarun = conepreprocess_wnm(datarun, varargin)
% CONEPREPROCESS_WNM    Setup datarun.stimulus for CONE_PREPROCESS on white noise movie
% usage: datarun = conepreprocess_wnm(datarun, opts)
%
% Sets up DATARUN.STIMULUS by loading the java_movie and then setting a 
% few additional fields.  OPTS are simply passed on to CONE_PREPROCESS.
%
% See also: CONE_PREPROCESS, CALC_CONE_INPUTS, GET_BINNED_SPIKES
%
% 2012-09 phli
%

if ~isfield(datarun, 'stimulus') || ~isfield(datarun.stimulus, 'java_movie')
    datarun = load_java_movie(datarun);
end
datarun.stimulus.refresh_time = datarun.stimulus.refresh_period/1000;
datarun.stimulus.stimobject = datarun.stimulus.java_movie;
datarun.stimulus.get_movie_frame = @get_movie_frame_from_java_movie;
datarun.stimulus.wc_preprocess = @wc_preprocess_javamovie_rawmovie;

if ~isfield(datarun, 'duration')
    datarun = load_neurons(datarun, 'load_spikes', []);
end
datarun.stimulus.duration = datarun.duration;

datarun = cone_preprocess(datarun, varargin{:});


function frame = get_movie_frame_from_java_movie(java_movie, i)
frame = double(java_movie.getFrame(i-1).getBuffer());
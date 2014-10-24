function datarun = load_java_movie(datarun, movie_xml_path, triggers)
% datarun = load_java_movie(datarun,movie_xml_path)
%
% arguments:     datarun - datarun struct with fields datarun.triggers, datarun.sampling_rate
%                              if these do not exist, then load_neurons is called
%         movie_xml_path - full path to xml file
%                           if not specified, used datarun.names.movie_xml_path
%
% outputs:     datarun - datarun struct with fields:
%
%           datarun.stimulus.java_movie      - java movie object
%           datarun.names.movie_xml_path     - path to xml movie
%
%  NOTE: if movie_xml_path does not specify a legitimate movie xml file, Java will cause matlab to crash
%
%
% 2009-09  gauthier
% 2012-07  phli, updated to load movie from .movie and .globals file if movie_xml_path is not given
% 2012-09  phli, can specify triggers instead of just using datarun.triggers
%

% Copy input path to datarun?
if nargin > 1 && ~isfield(datarun.names,'movie_xml_path')
    datarun.names.movie_xml_path = movie_xml_path;
end

% Use datarun's version?
if isfield(datarun.names, 'movie_xml_path') && nargin < 2
    movie_xml_path = datarun.names.movie_xml_path;
end


if exist('movie_xml_path', 'var')
    % Assume old kludgy method is desired...
    
    % ensure triggers and sampling rate are loaded
    if nargin < 3
        if ~isfield(datarun, 'triggers') || ~isfield(datarun, 'sampling_rate')
            datarun = load_neurons(datarun,'load_spikes,',[]);
        end
        triggers = datarun.triggers;
    end
    
    % show progress
    fprintf('\nComputing movie %s...',movie_xml_path);
    start_time = clock; % note when it started
    
    % ensure file exists, and ends in ".xml"
    if ~strcmpi(movie_xml_path(end-3:end),'.xml')
        error('movie xml file ''%s'' is not an xml file',movie_xml_path)
    end
    if ~exist(movie_xml_path,'file')
        error('movie xml file ''%s'' does not exist',movie_xml_path)
    end
    
    
    
    % generate movie
    movie = edu.ucsc.neurobiology.vision.matlab.Matlab.computeMovie(movie_xml_path, triggers*datarun.sampling_rate);
    
    % show duration
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
    
    % identify the stimulus in the java movie, and guess about parameters that were not specified
    the_stimulus = guess_stimulus(stimulus_from_java_movie(movie));
    stimsource = sprintf('XML movie %s',movie_xml_path);
    
elseif isfield(datarun.names, 'rrs_movie_path') && isfield(datarun.names, 'rrs_globals_path')
    % Newer preferred method
    movie = edu.ucsc.neurobiology.vision.matlab.Matlab.computeMovie(datarun.names.rrs_movie_path, datarun.names.rrs_globals_path);
    the_stimulus = stimulus_from_movie_and_globals(datarun);
    stimsource = '.movie and .globals files';
end

% load into datarun, and check for consistency 
datarun = sync_stimulus(datarun,the_stimulus,stimsource);
    
% incorporate into datarun
datarun.stimulus.java_movie = movie;


function stim = stimulus_from_movie_and_globals(datarun)
datarun = load_globals(datarun);
wnmovie = load_wnmovie(datarun);
stim = stimulus_from_globals(datarun.globals);
stim = setstructfields(stim, stimulus_from_wnmoviefile(wnmovie));
wnmovie.close;

function wnmovie = load_wnmovie(datarun)
READ = 0; % FIXME: Should use actual Java enum to get this value
wnmovie = edu.ucsc.neurobiology.vision.io.chunk.WhiteNoiseMovieFile(datarun.names.rrs_movie_path, READ);
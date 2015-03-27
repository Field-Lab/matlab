function [mvi extra] = load_movieAH(mdf_file, triggers)

% HACK for 20000 hz.  rather than puling param from Vision

% LOAD_MOVIE      Load a white noise movie
%
% usage : mvi = load_movie(xml_file, triggers);
%
% arguments:   xml_file - full path to movie xml file
%              triggers - trigger times from recording
%
% outputs:           mvi - Java object containing movie
%
%
% example: 
%    sta = load_movie(triggers, 'example.xml');
%
% shlens 2007-08-29
% tamachado 2009-06-15

% This AH version is just for the move form Salk to Stanford
% want the matlab portion of GLM code as self-contained as possible
% AH has not responsibililty in writing this!

% import java classes
import('edu.ucsc.neurobiology.vision.util.*');
import('edu.ucsc.neurobiology.vision.matlab.*');

% grab sampling frequency (in Hz)

%sampling_frequency = VisionParams.SAMPLES_PER_MILLISECOND * 1000;
sampling_frequency = 20000;

% check directory
if ~exist(mdf_file,'file')
  error(['generate_movie: movie xml file does not exist']);
end

% convert to integers
triggers = int32(triggers * sampling_frequency);

% generate buffered version of movie
mvi= Matlab.computeMovie(mdf_file, triggers);

% remove speedy buffer for simplicity
mvi = mvi.getUnderlyingMovie;

extra.height = double(mvi.getHeight);
extra.width = double(mvi.getWidth);    
extra.refreshtime = mvi.getRefreshTime; 
extra.nframes = mvi.size;
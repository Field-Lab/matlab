function [mvi extra] = load_movie(mdf_file, triggers)
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


% import java classes
import('edu.ucsc.neurobiology.vision.util.*');
import('edu.ucsc.neurobiology.vision.matlab.*');

% grab sampling frequency (in Hz)
sampling_frequency = VisionParams.SAMPLES_PER_MILLISECOND * 1000;

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
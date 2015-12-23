function [sta] = compute_sta_(spikes, triggers, vision_movie, varargin)
% COMPUTE_STA      Compute the spike triggered average
%
% usage : sta = compute_sta_(spikes, triggers, [vision_movie], <params>) 
%
% outputs:           sta - Java object containing entire STA.
%                          See vision.stimulus.STA
%
% optional params, their default values, and what they specify:
%
% depth           30               STA depth in frames
%
% offset          -1               How far past the spike the STA goes;
%                                  defaults to -1, which should have the
%                                  spike in the last frame.  0 would give 
%                                  an extra frame after the spike at the 
%                                  end.  (May not be correct for 
%                                  interval > 1.)
%
% See also COMPUTE_VISION_MOVIE

p = inputParser;

% specify list of optional parameters
p.addParamValue('depth', 30);
p.addParamValue('offset', -1);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% import java classes
import('edu.ucsc.neurobiology.vision.util.*');
import('edu.ucsc.neurobiology.vision.matlab.*');

% grab sampling frequency
sampling_frequency = VisionParams.samplesPerMillisecond * 1000;

% convert to integers
spikes   = int32(spikes   * sampling_frequency);
triggers = int32(triggers * sampling_frequency);

% calculate STA
sta = Matlab.computeSTA(spikes, triggers(1), vision_movie, params.depth, params.offset);
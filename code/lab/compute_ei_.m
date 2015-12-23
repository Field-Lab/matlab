function [ei position] = compute_ei_(spikes, data_path, varargin)
% COMPUTE_EI      Compute the electrophysiological image
%
% usage : [ei, position] = compute_ei(spikes, data_path, <params>)
%
% arguments:     spikes - 1xN list of spike times (s)
%             data_path - full directory path of raw data
%
%
% optional params, their default values, and what they specify:
%
% nlPoints           30          number of samples before each spike to use
% nrPoints           30          number of samples after each spike to use
% nSpikes            30          spikes to average over in EI computation
%
% example: [ei, position] = compute_ei(spikes, ...
%            '/Volumes/345/2005-04-26-1/data006/', 'nlPoints', 30);
%
% shlens 2006-07-23
% tamachado 2009-06-15
%greschner 2009-07-06

p = inputParser;

% specify list of optional parameters
p.addParamValue('nlPoints', 20);
p.addParamValue('nrPoints', 60);
p.addParamValue('nSpikes', length(spikes));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% import java classes
import('java.io.*');
import('edu.ucsc.neurobiology.vision.io.*');
import('edu.ucsc.neurobiology.vision.util.*');
import('edu.ucsc.neurobiology.vision.matlab.*');
import('edu.ucsc.neurobiology.vision.electrodemap.*');

% grab sampling frequency
sampling_frequency = VisionParams.samplesPerMillisecond * 1000;

% convert to integers
spikes = int32(spikes * sampling_frequency);

% check directory
if ~exist(data_path, 'dir')
  error(['compute_ei: bad directory path']);
end

% compute ei's
ei = Matlab.computeElectrophysiologicalImage(data_path, spikes, params.nlPoints, params.nrPoints, params.nSpikes);

% throw out the error
ei = squeeze(ei(1,:,:));

% throw out the triggers
ei = ei(2:end,:);

% grab paramters
electrodes = size(ei,1); 
frames = size(ei,2);

% grab array ID
file = File(data_path);
rawData = RawDataFile(file);
arrayID = int32(rawData.getHeader.getArrayID);
rawData.close();
clear rawData file;

% create java electrode map
electrodeMap = ElectrodeMapFactory.getElectrodeMap(arrayID);

% grab positions
for i=1:electrodes
  position(i,1) = electrodeMap.getXPosition(i);
  position(i,2) = electrodeMap.getYPosition(i);
end

if 0
    % zero disconnected electrodes
    for i=1:electrodes
      if electrodeMap.isDisconnected(i) == 1
        i
        ei(i,:) = 0;
      end
    end 
end

% cut off unused electrode
%ei = ei(1:electrodes-1,:);

% restore double
ei = double(ei);
position = double(position);












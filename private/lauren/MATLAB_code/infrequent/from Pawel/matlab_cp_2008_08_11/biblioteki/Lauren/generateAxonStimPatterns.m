function [electrodes Array] = generateAxonStimPatterns(centerElectrode)

% argument: a single electrode on the 61 electrode array
% returns: 
%
%   electrodes, an array of electrode numbers that are used in the
%   generated patterns
%
%   Array, an array of patterns, in which the first dimension corresponds
%   with an electrode in "electrodes" and the second dimension corresponds
%   with the pattern number.
%
%   see help for generate2ElectrodePatternsInClusterLockedAmp for details
%   about patterns generated

electrodes = getCluster(centerElectrode);

Array = generate2ElectrodePatternsInClusterLockedAmp(electrodes);
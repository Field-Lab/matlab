function spikeTimesInWindow = getNeuronSpikeTimes(pathToNeuronsFile, cellIDs, timeWindowStart, timeWindowEnd)

% retrieves spike times for a set of cells within a window of time from a neurons file, relative to
% the beginning of time window and in sample units
%
% arguments:
%   pathToNeuronsFiles, self-explanatory (should be a string)
%   cellIDs, a vector of cell IDs (determined from VISION analysis)
%   timeWindowStart, a scalar defining the beginning of the time window (in samples)
%   timeWindowEnd, a scalar defining the end of the time window (in samples)
%
% returns:
%   spikeTimesInWindow, a cell array of vectors denoting spike times of the chosen cells, as 
%   measured in samples from timeWindowStart
%
% author: Lauren Hruby
% created: 2008-08-15


%% for testing
% pathToNeuronsFile = '/Volumes/Lee/Analysis/Lauren/2008-04-22-2/data001-map/data001-map.neurons';
% cellIDs = [95 97 260 487];
% timeWindowStart = 2;
% timeWindowEnd = 3;

%%
spikeTimes = load_rrs_neurons(pathToNeuronsFile, cellIDs);

spikeTimesInWindow = cell(1, length(cellIDs));

for i = 1:length(cellIDs)
    firstSpikeIndex = find(spikeTimes{i}*20000 >= timeWindowStart, 1, 'first');
    lastSpikeIndex = find(spikeTimes{i}*20000 < timeWindowEnd, 1, 'last');
%     firstSpikeIndex = find(spikeTimes{i} >= timeWindowStart/20000, 1, 'first'); %time in seconds
%     lastSpikeIndex = find(spikeTimes{i} <= timeWindowEnd/20000, 1, 'last');
    spikeTimesInWindow{i} = 20000*(spikeTimes{i}(firstSpikeIndex:lastSpikeIndex) - timeWindowStart*ones(size(spikeTimes{i}(firstSpikeIndex:lastSpikeIndex))));
end
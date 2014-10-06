function [electrodes Array clusterIDs smallAmpPatterns] = generatePatternClusterStim2Secondaries512Wrapper(centerElectrodes, relAmps, varargin)


% usage:  [electrodes Array clusterIDs] = generatePatternClusterStim2Secondaries512Wrapper(centerElectrodes, relAmps)
%
% generates arrays that will be saved to create 'electrodes' and 'patterns' stimulus files for
% STIM512
%
% arguments:  
%           centerElectrodes - vector of electrodes that will be centers of clusters (clusters cannot be overlapping)
%                              ***clusters on edge of array are allowable
%             relAmps.normal - secondary electrode amplitudes, relative to the primary electrode
%                              amplitude, to be used in all combinations of 2 and 3 electrodes
%          relAmps.pairsOnly - secondary electrode amplitudes, relative to the primary
%                              electrode amplitude, to be used in combinations of 2 electrodes only (primary + 1 secondary)
%     relAmps.pairsOnlySmall - same as .pairsOnly, except that these patterns are played in a separate movie chunk from the rest to avoid
%                              rounding issues that occur when very low and very high currents are used on the same electrode in the same movie chunk
%
% outputs:        electrodes - a vector of electrode numbers that are used in the generated patterns
%                      Array - an array of patterns, in which the first dimension corresponds
%                              with an electrode in "electrodes" and the second dimension corresponds
%                              with the pattern number.
%                 clusterIDs - a vector of integers specifying which patterns correspond with
%                              which clusters; each cluster is given a unique integer value, and the vector
%                              indeces correspond with the pattern numbers (second dimension of Array)
%           smallAmpPatterns - binary vector that flags all patterns that
%                              need to be played in a separate movie chunk to avoid rounding errors
%
% author: Lauren Hruby Jepson, SNL-E
% last checked 2010-10-13
% added smallAmpThresh/smallAmpPatterns 2011-02-28
%
% 2011-04-01
% removed smallAmpThresh, replaced with relAmps.pairsOnlySmall, which
% specifies which secondary electrode amplitudes specifically to separate
% into different movie chunk
%
% Edited by Lauren Grosberg for STIM512 system, 2014-07-07

p = inputParser;

p.addRequired('centerElectrodes', @isnumeric)
p.addRequired('relAmps', @isstruct)

p.addParamValue('includeAnodalPrimaries', false, @islogical)

p.parse(centerElectrodes, relAmps, varargin{:})

includeAnodalPrimaries = p.Results.includeAnodalPrimaries;


% find nearest-neighboring electrodes (electrodes in each cluster)
clusterElectrodes = cell(1, length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodes{i} = getCluster512(centerElectrodes(i));
    if length(clusterElectrodes{i})<7
        h = warndlg('One of the chosen center electrodes is on the edge of the array.  Just FYI.');
        boxPos = get(h, 'position');
        set(h, 'position', [400 500 boxPos(3) boxPos(4)])
    end
end

%check if any clusters of electrodes overlap
if(clusterOverlapCheck512(clusterElectrodes))
    error('Clusters overlap.  Aborting.')
end

arrayForEachCluster = cell(1, length(centerElectrodes));
smallAmpPatternsEachCluster = cell(1, length(centerElectrodes));

for i = 1:length(centerElectrodes)
    if ~includeAnodalPrimaries
        [arrayForEachCluster{i} smallAmpPatternsEachCluster{i}] = generatePatternClusterStim2Secondaries(clusterElectrodes{i},...
            centerElectrodes(i), relAmps);
    else
        [arrayForEachCluster{i} smallAmpPatternsEachCluster{i}] = generatePatternClusterStim2SecondariesAnPrim512(clusterElectrodes{i},...
            centerElectrodes(i), relAmps);
    end
end

%concatenate patterns and electrodes and generate vector of cluster IDs
[Array electrodes clusterIDs smallAmpPatterns] = concatenatePatterns(arrayForEachCluster, clusterElectrodes, smallAmpPatternsEachCluster);



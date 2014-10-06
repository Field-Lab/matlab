function [catArray catElectrodes clusterIDs smallAmpPatterns] = concatenatePatterns(arrayForEachCluster, clusterElectrodes, smallAmpPatternsEachCluster)

% concatenates a set of arrays so that none of the columns or rows of one array are in the columns
% or rows of another, with zeros filling up extra space
%
% e.g. [1 2 3      and      [7 8       become     [1 2 3 0 0
%       4 5 6]               9 10]                 4 5 6 0 0
%                                                  0 0 0 7 8
%                                                  0 0 0 9 10]
%
% concatenates sets of electrodes in clusterElectrodes into a single vector of electrodes
%
% outputs:
%         catArray   concatenated array (specifies stimulus patterns; electrode x pattern number)
%    catElectrodes   concatenated vector of electrode numbers
%       clusterIDs   vector, with indices corresponding to second dimension of catArray, specifying
%                    a unique integer for all of the patterns associated with a given cluster
%    
%
%
% length of each clusterElectrodes vector should be equal to first dimension of corresponding
% arrayForEachCluster
%
% author: Lauren (SNL-E)
% updated 2009-08-20


%% for testing
% 
% arrayForEachCluster{1} = rand(5,6);
% arrayForEachCluster{2} = rand(10,7);
% arrayForEachCluster{3} = rand(2,4);
% 
% clusterElectrodes{1} = randperm(5);
% clusterElectrodes{2} = randperm(10);
% clusterElectrodes{3} = randperm(2);

%% new version

clusterIDs = [];

% measuring dimensions of each cluster's array
nClusters = length(arrayForEachCluster);
clusterDims = zeros(nClusters,2);
dimSum1 = 0; dimSum2 = 0;
for i = 1:nClusters
    clusterDims(i,:) = size(arrayForEachCluster{i});
    dimSum1 = dimSum1 + clusterDims(i,1); %dimension 1 of concatenated array
    dimSum2 = dimSum2 + clusterDims(i,2); %dimension 2 of concatenated array
    
    clusterIDs = [clusterIDs i*ones(1,clusterDims(i,2))]; %#ok<AGROW>
end

catArray = zeros(dimSum1, dimSum2);

% placing each cluster's array in appropriate location
dim1RunTotal = 1;
dim2RunTotal = 1;
for i = 1:nClusters
    catArray(dim1RunTotal : dim1RunTotal+clusterDims(i,1)-1, dim2RunTotal : dim2RunTotal+clusterDims(i,2)-1) = ...
        arrayForEachCluster{i};
    dim1RunTotal = dim1RunTotal+clusterDims(i,1);
    dim2RunTotal = dim2RunTotal+clusterDims(i,2);
end

catElectrodes = [clusterElectrodes{1}];
for i = 2:length(clusterElectrodes)
    catElectrodes = [catElectrodes clusterElectrodes{i}]; %#ok<AGROW>
end

smallAmpPatterns = smallAmpPatternsEachCluster{1};
for i = 2:length(smallAmpPatternsEachCluster)
    smallAmpPatterns = [smallAmpPatterns smallAmpPatternsEachCluster{i}]; %#ok<AGROW>
end
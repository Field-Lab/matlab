function n = add(n, e, cluster)
%ADD Add an electrode to the neurons file.
%   Provide output from clustering algorithm
%   and all the acceptable clusters present on electrode e will
%   be written to the neurons file.
%
%   tamachado@salk.edu 2/4/08

data = toStruct(n.data{1}, e);

clustersAdded = zeros(1, length(cluster.assignments));

% Since assignments is a cell array, this will tell us how many 
% solutions we are currently saving out to neurons files
nSolutions = length(cluster.assignments);

for k = 1:nSolutions
    
    nClusters = cluster.nClusters(k);
    
    for c = 1:nClusters
        spikeTimes = data.spikeTimes(find(cluster.assignments{k} == c));
        nSpikes    = length(spikeTimes);

        if (nSpikes > n.minSpikes{1})
            n.nObject{k}.addNeuron(e, clustersAdded(k), spikeTimes, nSpikes);
            clustersAdded(k) = clustersAdded(k) + 1;
        end
    end
    
    disp(sprintf('%d clusters added to electrode %d.\n', clustersAdded(k), e))
end



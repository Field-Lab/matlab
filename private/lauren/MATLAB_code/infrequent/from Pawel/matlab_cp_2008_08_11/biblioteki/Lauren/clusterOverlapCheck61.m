function trueOrFalse = clusterOverlapCheck61(clusterElectrodes)

% checks to make sure clusters don't overlap
% argument should be a cell array containing an array of electrode numbers for each cluster
% returns 1 if there is overlap and 0 otherwise

trueOrFalse = 0;
for i = 1:64
    for j = 1:length(clusterElectrodes)-1
        for k = j+1:length(clusterElectrodes)
            if ~isempty(find(i==clusterElectrodes{j},1)) && ~isempty(find(i==clusterElectrodes{k},1))
                trueOrFalse = 1;
            end
        end
    end
end

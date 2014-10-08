function combs = nearestNeighborsFinder512(electrodes)

% arguments:
%   electrodes -- a vector of electrode numbers from a 512-electrode array
%
% given a set of a electrodes on the 512 array, this function determines
% which pairs of electrodes within the cluster are nearest-neighbors (30 or
% 60 microns apart for the 30 and 60 micron arrays, respectively), 
%
% output: an array of electrode pairs ("combinations"), with the first
% index representing the combination number, and the second index denoting
% the indeces (in "electrodes") of each of the two electrodes in each combination.

[xCoords, yCoords] = getElectrodeCoords512();

nElectrodes = length(electrodes);

combs = [];
%combinations = zeros(nCombinations, 2);

%combinationCount = 0;
for i = 1:nElectrodes
    for j = i+1:nElectrodes
        disp(['Electrode ' num2str(electrodes(i)) ' and ' num2str(electrodes(j)) ' are separated by '...
            num2str(norm([xCoords(electrodes(i)) - xCoords(electrodes(j)), yCoords(electrodes(i)) - yCoords(electrodes(j))]))]);
       if norm([xCoords(electrodes(i)) - xCoords(electrodes(j)), yCoords(electrodes(i)) - yCoords(electrodes(j))]) < 68
           combs = [combs; i j]; %#ok<AGROW>
           %combinationCount = combinationCount + 1;
           %combinations(combinationCount, 1) = i;
           %combinations(combinationCount, 2) = j;
       end
    end
end


% if combinationCount ~= nCombinations
%     error('The expected number of nearest-neighbor pairs does not match the number found by nearestNeighborsFinder')
% end

function combinations = nearestNeighborsFinder(electrodes, nCombinations)

% arguments:
%   electrodes -- a vector of electrode numbers from a 61-electrode array
%   nCombinations -- the expected number of nearest-neighbor pairs within
%       the set of electrodes (produces an error if this number is 
%       incorrect)
%
% given a set of a electrodes on the 61 array, this function determines
% which pairs of electrodes within the cluster are nearest-neighbors (30 or
% 60 microns apart for the 30 and 60 micron arrays, respectively), and
% returns an array of electrode pairs ("combinations"), with the first
% index representing the combination number, and the second index denoting
% each of the two electrodes in each combination.  The values in the array
% are the indeces of the paired electrodes in "electrodes."

[xCoords yCoords] = getElectrodeCoords61();

nElectrodes = length(electrodes);

combinations = zeros(nCombinations, 2);

combinationCount = 0;
for i = 1:nElectrodes
    for j = i+1:nElectrodes
       if norm([xCoords(electrodes(i)) - xCoords(electrodes(j)), yCoords(electrodes(i)) - yCoords(electrodes(j))]) < 2.1
           combinationCount = combinationCount + 1;
           combinations(combinationCount, 1) = i;
           combinations(combinationCount, 2) = j;
       end
    end
end

if combinationCount ~= nCombinations
    error('The expected number of nearest-neighbor pairs does not match the number found by nearestNeighborsFinder')
end

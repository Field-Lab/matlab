function [Times Patterns] = combineNeuronSpikeTimes(spikeTimes, patternNumbers)

nCells = length(spikeTimes);

nTimes = 0;
for i=1:nCells
    nTimes = nTimes + length(spikeTimes{i});
end
    
Times = zeros(1, nTimes);
Patterns = zeros(1, nTimes);

nextTimes = zeros(1, nCells); %next (unadded) spike time for each cell
nextIndeces = ones(1, nCells); %indeces of nextTimes

for i = 1:nTimes
    for j = 1:nCells
        if nextIndeces(j) > length(spikeTimes{j})
            nextTimes(j) = inf;
        else
            nextTimes(j) = spikeTimes{j}(nextIndeces(j));
        end
    end
    nextCellIndex = find(nextTimes == min(nextTimes));
    Patterns(i) = patternNumbers(nextCellIndex);
    Times(i) = spikeTimes{nextCellIndex}(nextIndeces(nextCellIndex));
    nextIndeces(nextCellIndex) = nextIndeces(nextCellIndex) + 1;
end
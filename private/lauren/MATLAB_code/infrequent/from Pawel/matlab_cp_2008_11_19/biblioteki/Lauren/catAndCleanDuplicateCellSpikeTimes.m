function spikeTimesClean = catAndCleanDuplicateCellSpikeTimes(pathToNeuronsFile, cells, timeWindow)

movingBarStartTime = timeWindow(1);
movingBarEndTime = timeWindow(2);

spikeTimesClean = cell(1,length(cells));

for i = 1:length(cells)
    spikeTimes = getNeuronSpikeTimes(pathToNeuronsFile, cells{i}, movingBarStartTime, movingBarEndTime);
    spikeTimesTotal = [];
    for j = 1:length(spikeTimes)
        spikeTimesTotal = [spikeTimesTotal; spikeTimes{j}];
    end
    spikeTimesTotal = sort(spikeTimesTotal);
    for j = length(spikeTimesTotal)-1:-1:1;
        if (spikeTimesTotal(j+1) - spikeTimesTotal(j) <= 10); %if two spike times are within 0.5 ms from eachother
            spikeTimesTotal(j+1) = []; % get rid of the second spike time
        end
    end
    spikeTimesClean{i} = spikeTimesTotal;
end


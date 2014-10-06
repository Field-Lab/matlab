function [successRates stimAmps] = getSuccessRatesFromLatencies(dataPath, patternNumbers, movieNumbers, latencies)


successRates = zeros(size(latencies));
stimAmps = zeros(size(latencies));

for i = 1:length(patternNumbers)
    for j = 1:length(movieNumbers)
        successRates(i,j) = sum(latencies{i,j}~=0)/length(latencies{i,j});
        stimAmps(i,j) = max(abs(getStimAmps(dataPath, patternNumbers(i), movieNumbers(j))));
    end
end


for j = 1:length(movieNumbers)
    stimAmps(1,j) = max(abs(getStimAmps(dataPath, patternNumbers(1), movieNumbers(j))));
end
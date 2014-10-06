function bootstrapNDrawsCompare(pathToData, patternNos, neuronID, varargin)

% a wrapper function to determine the effect of number of pulses on threshold standard deviations

p = inputParser;

p.addRequired('pathToData', @ischar)
p.addRequired('patternNos', @isnumeric)
p.addRequired('neuronID', @isnumeric)

p.addParamValue('bootstrapReps', 100, @isnumeric)
p.addParamValue('recalcAll', 0, @(x)any(x==[0 1]))

p.parse(pathToData, patternNos, neuronID, varargin{:})

nBootstrapReps = p.Results.bootstrapReps;
recalcAll = p.Results.recalcAll;



cd(pathToData)

stdevAllPulses = zeros(length(patternNos), 1);
meanAllPulses = zeros(length(patternNos), 1);
stdevNDraws = zeros(length(patternNos), 1);
meanNDraws = zeros(length(patternNos), 1);
actualThresh = zeros(length(patternNos), 1);

for i = 1:length(patternNos)
    disp(num2str(i))
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
    elecResp = temp.elecResp;


    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll, 'keepLogBased', 0);
    save(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))], 'elecResp')
    
    actualThresh(i) = elecResp.analysis.threshold;
    
    result = bootstrapThresh(elecResp, nBootstrapReps, 'nDraws', 25);
    stdevAllPulses(i) = result.stdevAllPulses;
    meanAllPulses(i) = result.meanAllPulses;
    stdevNDraws(i) = result.stdevNDraws;
    meanNDraws(i) = result.meanNDraws;
end


%%plotting results

figure
hold on
plot(stdevAllPulses./actualThresh, 'k')
plot(stdevNDraws./actualThresh, 'b')
plot(abs(meanAllPulses-meanNDraws)./actualThresh, 'r')
hold off
legend('standard deviation for full 50 pulses', 'standard deviation using first 25 pulses', 'abs(difference) between mean threshold using 50 pulses and mean using 25 pulses')
ylabel('normalized to actual threshold')
xlabel('trial')

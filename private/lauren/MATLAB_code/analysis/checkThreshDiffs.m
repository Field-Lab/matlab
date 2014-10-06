
neuronID = 91; 

primThresh = -1.2612;

patternNos = 1:36;
nPatterns = length(patternNos);

threshDiffs = zeros(36, 1);
for i = 1:nPatterns
    disp(num2str(i))
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
    elecResp = temp.elecResp;
    
    threshDiffs(i) = elecResp.analysis.constrainedSlopeLogThresh - primThresh;
end
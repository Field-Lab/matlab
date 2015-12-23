function [nPosSpikes, nNegSpikes, varSpikesElicitedPos, varSpikesElicitedNeg] = ...
    computeNSpikesElicitedInPSTH(PSTH, startTimeBaseline, endTimeBaseline)
% [nPosSpikes, nNegSpikes, varSpikesElicitedPos, varSpikesElicitedNeg] = ...
%     computeNSpikesElicitedInPSTH(PSTH, startTimeBaseline, endTimeBaseline)
% 
% Parameters:
%   - PSTH: a psth struct. Bin size is assumed to be 5ms. Fields are data, 
%   variances.
%   - startTimeBaseline: start time of the quiet period.
%   - endTimeBaseline: end time of the quiet period.
%
% Returns:
%   - nPosSpikes: number of spikes elicited
%   - nNegSpikes number of spikes suppressed
%   - varSpikesElicited: variance in the number of spikes elicited.
%

pThreshold = 0.01;

% Getting the baseline statistics.
[nSpikesBackground, semBackground, nTrials, nBins] = getNumberOfSpikes(...
    startTimeBaseline, endTimeBaseline, PSTH);
nSpikesExpected = nSpikesBackground/nBins; % Number of spikes expected per bin

% Estimate from SEM is biased, hence the correcting factor.
% Other assumption made is that what happens in a bin is
% independent from what happens in others, which is an
% approximation but should be close to true.
varSpikesExpected = nTrials/(nTrials-1)*semBackground*nTrials/nBins;

% Going over all the bins and finding the statistically significant
% ones, adding their contribution to the number of spikes elicited.
nPosSpikes = 0;
nNegSpikes = 0;
varSpikesElicitedPos = 0;
varSpikesElicitedNeg = 0;
for mm=3:length(PSTH.data)-1 % First bins are usually blanked or artifacts
    nSpikesBin = PSTH.data(mm);
    varBin = PSTH.variances(mm);

    tValue = abs((nSpikesBin-nSpikesExpected)/...
        sqrt(varSpikesExpected/nTrials + varBin/nTrials));
    pValue = 1-tcdf(tValue, nTrials-1);

    if pValue < pThreshold
        if nSpikesBin>nSpikesExpected % Stimulation
            nPosSpikes = nPosSpikes + (nSpikesBin - nSpikesExpected);
            varSpikesElicitedPos = varSpikesElicitedPos + varBin;
        else % Suppression
            nNegSpikes = nNegSpikes + (nSpikesBin - nSpikesExpected);
            varSpikesElicitedNeg = varSpikesElicitedNeg + varBin;
        end
        
    end
end

end % computeNSpikesElicitedInPSTH

function [nSpikes, sem, nTrials, nBins] = getNumberOfSpikes(startTime,endTime,PSTH)
% startTime, endTime specified in ms. 
% Assuming bin size of 5 ms in the PSTH, Fs = 20,000 Hz
% Returns the standard error of the mean.

x = (floor(startTime/5)+1):(floor(endTime/5));
nSpikes = sum(PSTH.data(x));
sem = sum(PSTH.variances(x))/PSTH.numberOfPulses; % Square of SEM. Addition because Poisson distr. assumption
nTrials = PSTH.numberOfPulses;
nBins = length(x);

end % getNumberOfSpikes
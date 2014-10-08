function [electrodes Array] = generateSpatialPatternStimPatterns(electrodes, stimAmpsFull, stimAmpsThresh)

% generates arrays that will be saved to create 'electrodes' and 'patterns' stimulus files for
% STIM64
%
% arguments
%    electrodes: vector of electrodes to be used in spatial pattern (each targets a different cell)
%
% returns
%   electrodes, an array of electrode numbers that are used in the
%   generated patterns (same as argument)
%
%   Array, an array of patterns, in which the first dimension corresponds
%   with an electrode in "electrodes" and the second dimension corresponds
%   with the pattern number.
%
%   ***all amplitudes specified in Array are relative to the maximum amplitude in stimAmpsFull
%
%
%   pattern 1: all electrodes stimulated simultaneously at amplitudes that give 100% response rate
%   patterns 2-(1+nElec): each electrode stimulated individually at amplitude that gives 100% response
%   rate
%
%   patterns (2+nElec)-(1+6*nElec): all electrodes stimulated simultaneously, with all but one of the
%   electrodes stimulated at amplitudes that give 100% reponses, and last electrode stimulated at
%   amplitude at one of 5 amplitudes around threshold
%
%   patterns (2+6*nElec)-(1+11*nElec): each electrode stimulated individually at amplitude near
%   threshold
%

nElec = length(electrodes);
maxAmp = max(stimAmpsFull);

Array = zeros(nElec, 11*nElec+1);
relativeAmps = [0.81 0.9 1 1.1 1.21]; %amplitudes, relative to threshold, to be tested to detect changes in response curve

for i = 1:nElec
    %all electrodes stimulated simultaneously at amplitudes that give 100% response rate
    Array(i, 1) = stimAmpsFull(i)/maxAmp;
    %each electrode stimulated individually at amplitude that gives 100% response rate
    Array(i, 1+i) = stimAmpsFull(i)/maxAmp;
    
    for j = 1:nElec
        for k = 1:5
            %all electrodes stimulated simultaneously at amplitudes that give 100% response rate, except for
            %one electrode which is stimulated at 5 different amplitudes around threshold
            if i==j
                Array(i, 1+nElec+5*(j-1)+k) = stimAmpsThresh(i)*relativeAmps(k)/maxAmp;
            else
                Array(i, 1+nElec+5*(j-1)+k) = stimAmpsFull(i)/maxAmp;
            end
            
            % each electrode stimulated alone near thresh
            if i==j
                Array(i, 1+6*nElec+5*(j-1)+k) = stimAmpsThresh(i)*relativeAmps(k)/maxAmp;
            end
        end
    end
end


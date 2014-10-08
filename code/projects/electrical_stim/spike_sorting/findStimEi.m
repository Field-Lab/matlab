function ei = findStimEi(dataTraces, elecResp, movieNumber, tempIndex, eiData, channelsToUse, artFromSubset)

% usage: ei = findStimEi(dataTraces, elecResp, movieNumber, tempIndex, eiData, channelsToUse)
%
% purpose: isolates portions of data representing spike of specified neuron (removes artifact and spikes of
% other neurons specified in elecResp.analysis.latencies).  Returns data vectors on chosen
% electrodes, positioning the spike minima at sample 11 of 26
%
%
% arguments:
%    dataTraces: nTraces x nElectrodes x nSamples
%    elecResp: struct containing information about data/analysis (see createElecRespStruct.m)
%    movieNumber: self-explanatory (actual movie number, not index)
%    tempIndex: index corresponding to template within elecResp.cells.active 
%       **if tempIndex == 0, signifies elecResp.cells.main (main neuron)
%    eiData: cell array (indeces corresponding to elecResp.cells.active+1) of arrays (channels,
%    samples), including all channels other than 0 (the TTL channel)
%    channelsToUse: vector of channels number to include in returned ei
%       **if there are no traces that are artifact-only (no identified spikes), findStimEi will
%       replace channelsToUse with the channels listed in elecResp.cells.goodElecs
%    artFromSubset: logical specifying whether to include (false) or
%    disclude (true) traces that are marked with an analysis flag when
%    calculating the artifact estimate
%
%
% returns:
%    ei = traces x channels x samples
%
% author: Lauren Hruby (SNL-E), early 2009
% last updated: 2009-06-16


movieIndex = find(elecResp.stimInfo.movieNos == movieNumber);

centerChannel = elecResp.cells.recElec;
goodChannels = elecResp.cells.goodElecs;
neuronIDs = [elecResp.cells.main elecResp.cells.active{movieIndex}];
if isfield(elecResp.analysis, 'estArtifact')
    estArtifact = elecResp.analysis.estArtifact{movieIndex};
end
    
if elecResp.stimInfo.nPulses(movieIndex) == 0 %value hasn't been set yet
    dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, movieNumber);
    elecResp.stimInfo.nPulses(movieIndex) = size(dataTraces, 1);
end
nTraces = elecResp.stimInfo.nPulses(movieIndex);

latencies = [elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];

nTemplates = length(neuronIDs);

failuresBin = (latencies == 0);
failuresBin = min(failuresBin, [], 2);

tempIndex = tempIndex + 1; %because latencies and neuronIDs include active neurons concatenated to main neuron

if artFromSubset %when calculating artifact, only use unflagged failure traces
    unflaggedBin = elecResp.analysis.details.analysisFlags{movieIndex}==0;
else
    unflaggedBin = true(size(elecResp.analysis.details.analysisFlags{movieIndex}));
end

failuresBin = logical(failuresBin.*min(unflaggedBin,[],2));


if max(failuresBin) %some traces have no identified spikes
    estArtifact = reshape(mean(dataTraces(failuresBin, channelsToUse, :), 1), length(channelsToUse), []);
elseif exist('estArtifact','var')
    channelsToUse = goodChannels;
else
    ei = [];
    return
end
nChannels = length(channelsToUse);

ei = zeros(nTraces, nChannels, 26); 

for k = 1:nTraces
    if latencies(k, tempIndex)

        % subtracts ei of other neurons that have been found in same trace
        for i = 1:nTemplates
            if i~=tempIndex && latencies(k,i)~=0 %if there is a spike from a different neuron
                otherEi{1} = eiData{i}(channelsToUse, :); %#ok<AGROW>

                otherEiMinPos = find(eiData{i}(centerChannel,:) == min(eiData{i}(centerChannel,:)));
                offset = latencies(k, i) - otherEiMinPos;
                trace = reshape(dataTraces(k, channelsToUse, :), length(channelsToUse), []);
                dataTraces(k, channelsToUse, :) = subtractWithShifts(trace, otherEi, offset);
            end
        end

        latency = latencies(k, tempIndex); %location of minimum in trace

        subtractedTrace = reshape(dataTraces(k, channelsToUse, :), length(channelsToUse), []) - estArtifact;

        if mod(latency,1) %need to interpolate
            nSamples = size(subtractedTrace, 2);
            for m = 1:nChannels
                subtractedTrace(m, :) = pchip(1:nSamples, squeeze(subtractedTrace(m,:)),...
                    (1:nSamples) + mod(latency,1));
            end
            latency = floor(latency);
        end
        if latency > 10 && size(subtractedTrace,2) >= latency+15
            ei(k, :, :) = subtractedTrace(:, latency-10 : latency+15);
        elseif latency > 10
            missingSamp = latency+15 - size(subtractedTrace,2);
            ei(k, :, 1:end-missingSamp) = subtractedTrace(:, latency-10 : end);
        else
            missingSamp = 11 - latency;
            ei(k, :, 1+missingSamp : 26) = subtractedTrace(:, 1 : 26-missingSamp);
        end
    end
end

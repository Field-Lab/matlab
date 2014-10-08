function elecResp = templateMatchClustering1Trace(elecResp, movieNumber, traceNo, flags)

% shorter version of templateMatchClustering, designed to just reanalyze a single trace
%
% differences from elecResp:
%   model artifact is the estimated artifact from previous analysis
%   
%
%

%% argument parsing

p = inputParser;

p.addRequired('elecResp', @isstruct)
p.addRequired('movieNumber', @isnumeric)
p.addRequired('tracNo', @isnumeric)
p.addRequired('flags', @isnumeric)

p.parse(elecResp, movieNumber, traceNo, flags)

movieIndex = find(movieNumber == elecResp.stimInfo.movieNos);

shiftStart = elecResp.analysis.details.tempOffsetWindow{movieIndex}(1);
shiftEnd = elecResp.analysis.details.tempOffsetWindow{movieIndex}(2);
shiftStep = elecResp.analysis.details.tempOffsetStep{movieIndex};
nShifts = (shiftEnd - shiftStart)/shiftStep + 1;

startFit = elecResp.analysis.details.residCalcWindow{movieIndex}(1); 
endFit = elecResp.analysis.details.residCalcWindow{movieIndex}(2);

% extracting info from elecResp
dataPath = elecResp.names.data_path;
patternNumber = elecResp.stimInfo.patternNo;
pathToEi = elecResp.names.rrs_ei_path;
neuronID = [elecResp.cells.main elecResp.cells.active{movieIndex}];
centerChannel = elecResp.cells.recElec;
goodChannels = elecResp.cells.goodElecs;
prevLatencies = [elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];
estArtifact = elecResp.analysis.estArtifact{movieIndex};

centerChannelIndex = find(goodChannels == centerChannel);
nElecs = length(goodChannels);

%%

% if any of the analysis flags is 1 (skip reanalysis), checks to make sure data has been previously
% analyzed
if any(flags==1) && isempty(prevLatencies)
    warnH = warndlg(['Data must be previously analyzed (elecResp.analysis.latencies must not be empty)'...
        'in order to skip the reanalysis of some traces/templates']);
    uiwait(warnH)
    return
end

%% loads data

%returns data traces: nTraces x nElectrodes x nSamples
dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, patternNumber, movieNumber, 99999);

nTemplates = length(neuronID);

templates = cell(nTemplates, 1);
templateMinPos = zeros(nTemplates, 1);
% if isfield(elecResp.cells, 'allEIs')
    for i = 1:nTemplates
        templates{i} = elecResp.cells.allEIs{elecResp.cells.all == neuronID(i)}(goodChannels, :);
        templateMinPos(i) = find(squeeze(templates{i}(centerChannelIndex,:)) ==...
            min(squeeze(templates{i}(centerChannelIndex,:))));
    end
% else
%     eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
%     %ei = cell(nTemplates, 1);
%     for i = 1:length(neuronID)
%         ei = eiFile.getImage(neuronID(i)); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
%         templates{i} = reshape(ei(1, goodChannels + 1, :), length(goodChannels), []);
%         
%         templateMinPos(i) = find(squeeze(templates{i}(centerChannelIndex,:)) ==...
%             min(squeeze(templates{i}(centerChannelIndex,:))));
%     end
%     clear eiFile
% end


%%
if any(flags == 1)
    offsets = shiftStart:shiftStep:shiftEnd;
    prevOffsetInd = zeros(1, nTemplates);
    for i = 1:nTemplates
        prevOptOffset = prevLatencies(i) - templateMinPos(i);
        prevOffsetInd(i) = find(offsets == prevOptOffset);
    end
else
    prevOffsetInd = [];
end



%% 

optOffsets = zeros(1, nTemplates);

%%% subtracts template waveform from data at different latencies
testTrace = reshape(dataTraces(traceNo, goodChannels, :), length(goodChannels), []);
[subtractedTraces offsets] = subtractWithShifts(testTrace, templates, shiftStart, shiftEnd, shiftStep);

%%% finding error between model artifact and traces after subtraction
errors = calcErrors(subtractedTraces, 'hypArtifact', startFit, endFit, estArtifact, flags, prevOffsetInd);

optErrors = min(reshape(errors,[],1));
minInd = find(errors == optErrors,1);

% a horrible kluge to make ind2sub actually WORK
if nTemplates < 2
    minSub(1) = find(errors == optErrors, 1);
else
    if nTemplates < 3
        [minSub(1) minSub(2)] = ind2sub(size(errors), find(errors == optErrors, 1));
    else
        if nTemplates < 4
            [minSub(1) minSub(2) minSub(3)] = ind2sub(size(errors), find(errors == optErrors, 1));
        else
            [minSub(1) minSub(2) minSub(3) minSub(4)] = ind2sub(size(errors), find(errors == optErrors, 1));
        end
    end
end

% removing error minima that lie at smallest and largest offset values
noMoreExtremeMins = 0;

while ~noMoreExtremeMins
    if any(minSub==2) %if lowest offset gives lowest error, set it and the offsets within the next 10 samples of that template to infinity
        errors(minInd) = inf;
        badDim = find(minSub==2, 1);
        for i = 1:(10/shiftStep)
            afterMinInd = minInd + i*size(errors,1)^(badDim-1);
            errors(afterMinInd) = inf;
        end
    elseif any(minSub==(nShifts+1)) %if highest offset gives lowest error, set it and the offsets within the previous 10 samples of that template to infinity
        errors(minInd) = inf;
        badDim = find(minSub == nShifts+1, 1);
        for i = 1:(10/shiftStep)
            beforeMinInd = minInd - i*size(errors,1)^(badDim-1);
            errors(beforeMinInd) = inf;
        end
    elseif all(errors==inf) %no local minima could be found at required distance from edges
        errordlg('no local minima could be found further than 10 samples from boundaries');
        return
    else
        noMoreExtremeMins = 1;
    end

    optErrors = min(reshape(errors,[],1));
    minInd = find(errors == optErrors,1);

    % a horrible kluge to make ind2sub actually WORK
    if nTemplates < 2
        minSub(1) = minInd;
    else
        if nTemplates < 3
            [minSub(1) minSub(2)] = ind2sub(size(errors), minInd);
        else
            if nTemplates < 4
                [minSub(1) minSub(2) minSub(3)] = ind2sub(size(errors), minInd);
            else
                [minSub(1) minSub(2) minSub(3) minSub(4)] = ind2sub(size(errors), minInd);
            end
        end
    end
end

% getting offset values corresponding to error mimimum
for i = 1:nTemplates
    optOffsets(i) = offsets(minSub(i));
end


%% determines latencies

latencies = zeros(1, nTemplates);

if all(isnan(optOffsets))
else
    for j = 1:nTemplates
        if ~isnan(optOffsets(j))
            latencies(j) = optOffsets(j) + templateMinPos(j);
        end
    end
end

%% calculates artifact based on final classification

nSuccesses = 0;
actualArtifact = zeros(nElecs, size(dataTraces,3)); %average of failure traces, based on full algorithm
latenciesFull = prevLatencies;
latenciesFull(traceNo, :) = latencies;

failures = 0;
for i = 1:size(dataTraces, 1)
    if all(latenciesFull(i,:) == 0)
        actualArtifact = actualArtifact + reshape(dataTraces(i, goodChannels, :), length(goodChannels), []);
        failures = failures + 1;
    elseif latenciesFull(i, 1)
        nSuccesses = nSuccesses + 1;
    end
end

actualArtifact = actualArtifact/failures;
successRate = nSuccesses/size(dataTraces,1);


%% updating elecResp

elecResp.analysis.latencies{movieIndex} = latenciesFull(:,1);
if ~isempty(elecResp.cells.active{movieIndex})
    elecResp.analysis.otherLatencies{movieIndex} = latenciesFull(:,2:end);
end
elecResp.analysis.erfCurrent = 0;
elecResp.analysis.finalized(movieIndex) = 0;
elecResp.analysis.successRates(movieIndex) = successRate;

elecResp.analysis.details.analysisFlags{movieIndex}(traceNo, :) = flags;

if ~any(any(isnan(actualArtifact))) && max(max(abs(actualArtifact)))
    elecResp.analysis.estArtifact{movieIndex} = actualArtifact;
end

end






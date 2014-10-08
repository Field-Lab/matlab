function [latencies latenciesInitial actualArtifact] = spikeSubtractionLinkage(dataPath, patternNumber, movieNumber,...
    pathToEi, neuronID, centerChannel)


%% additional parameters

linkageWindow = [10 40];
linkThresh = 30;
nBranches = 2;

startLinFit = 10;
endLinFit = 40;


shiftStart = -10; %(spike starts at about sample 20 in template)
shiftEnd = 10;
shiftStep = 0.25;


saveFigures = 1;

channelRadius = 1;


%%


%returns data traces: nTraces x nElectrodes x nSamples
dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, patternNumber, movieNumber);

nTraces  = size(dataTraces, 1);
nSamples = size(dataTraces, 3);

% final output
latencies = zeros(1, nTraces); %if no spike, = 0, otherwise, distance from pulse start to minimum of spike

%% gets 'template' ei from .ei file

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples

clear eiFile

%% attempts to find artifact based on linkage analysis

linkData = dataTraces(:, centerChannel, :);
templates = squeeze(ei(1, centerChannel + 1, :));

estArtifact = estArtWithLinkage(linkData, linkageWindow, linkThresh, nBranches, templates);


%% loops through traces in data

optOffsets = zeros(1, size(dataTraces,1));
optSubtracted = cell(1, size(dataTraces,1));
optErrors = zeros(1, size(dataTraces,1));
allErrors = cell(1, size(dataTraces,1));

errorsNoSub = zeros(1,size(dataTraces,1));

template = squeeze(ei(1, centerChannel + 1, :));
templateMinPos = find(template == min(template));

for k = 1:size(dataTraces,1)
    
    disp(['computing for trace ' num2str(k)])

    %% subtracts template waveform from data at different latencies
    testTrace = squeeze(dataTraces(k, centerChannel, :));
    [subtractedTraces offsets] = subtractWithShifts(testTrace, template, shiftStart, shiftEnd, shiftStep);
    
    %% fitting exponential to traces before subtraction
    if artifactSubtraction && ~expFitOverride
        errorsNoSub(k) = norm(testTrace(startLinFit:endLinFit) - mean(testTrace(startLinFit:endLinFit)));
    elseif usePrevArtifact == 1 && length(previousArtifact) > 1 && ~isnan(max(previousArtifact))
        %errorsNoSub(k) = artFitter(testTrace(startArtFit:endArtFit), previousArtifact, [startArtFit, endArtFit]);
        errorsNoSub(k) = norm(testTrace(startArtFit:endArtFit) - previousArtifact(startArtFit:endArtFit));
    elseif useSumOfExp
        errorsNoSub(k) = expFitter2(testTrace(startExpFit:endExpFit));
    else
        errorsNoSub(k) = expFitter(testTrace(startExpFit:endExpFit));
    end


    %% fitting exponential to the artifact (traces after subtraction) 

    %params = zeros(3, length(subtractedTraces));
    projections = cell(1, length(subtractedTraces));
    errors = zeros(1, length(subtractedTraces));
    h.axes = cell(1, length(subtractedTraces));
    
    for i = 1:length(subtractedTraces)
        if artifactSubtraction && ~expFitOverride
            errors(i) = norm(subtractedTraces{i}(startLinFit:endLinFit)...
                - mean(subtractedTraces{i}(startLinFit:endLinFit)));
        elseif usePrevArtifact == 1 && length(previousArtifact) > 1 && ~isnan(max(previousArtifact))
            errors(i) = norm(subtractedTraces{i}(startArtFit:endArtFit) - previousArtifact(startArtFit:endArtFit));
        elseif useSumOfExp
            [errors(i) projections{i}] = expFitter2(subtractedTraces{i}(startExpFit:endExpFit));
        else
            errors(i) = expFitter(subtractedTraces{i}(startExpFit:endExpFit));
        end

    end
    
    optOffsets(k) = offsets(find(errors == min(errors), 1));
    optSubtracted{k} = subtractedTraces{find(errors==min(errors), 1)};
    optErrors(k) = min(errors);
    allErrors{k} = errors;

    if generatePlots
        plotWithSlider(h.axes)

        figure
        plot(offsets, errors)

        title(['minimum occurs at offset = ' num2str(optOOffsets(k))])
    end
end


%% chooses error that is lower from optErrors and errorsNoSub

lowerExpError = zeros(nTraces, 1);
lowerExpErrorTrace = cell(nTraces, 1);
successesExpFit = zeros(nTraces, 1);
actualArtifactExpOnly = zeros(nSamples, 1);
latenciesInitial = zeros(nTraces, 1);

for i = 1:nTraces
    if optErrors(i) <= errorsNoSub(i) && optOffsets(i) ~= shiftStart && optOffsets(i) ~= shiftEnd
        lowerExpError(i) = optErrors(i);
        lowerExpErrorTrace{i} = optSubtracted{i};
        successesExpFit(i) = 1;
        latenciesInitial(i) = optOffsets(i) + templateMinPos;
    else
        lowerExpError(i) = errorsNoSub(i);
        lowerExpErrorTrace{i} = squeeze(dataTraces(i, centerChannel, :));
        actualArtifactExpOnly = actualArtifactExpOnly + squeeze(dataTraces(i, centerChannel, :));
    end
end
actualArtifactExpOnly = actualArtifactExpOnly/sum(~successesExpFit);
successesExpFit = cast(successesExpFit, 'logical'); %not sure why this is necessary, but get an error otherwise

%% creates hypothetical artifact from best fits to exponential (from subtracted and non-subtracted)

meanExpError = mean(lowerExpError);
artifactSum = zeros(length(optSubtracted{1}), 1);
nArtifacts = 0;

for i = 1:nTraces
    if lowerExpError(i) < meanExpError
        artifactSum = artifactSum + lowerExpErrorTrace{i};
        nArtifacts = nArtifacts + 1;
    end
end

artifactMean = artifactSum / nArtifacts;

%% redoing error minimization of template subtraction at different latencies based on hyp. artifact

optOffsetsFull = zeros(nTraces, 1);
optSubtractedFull = cell(nTraces, 1);
optErrorsFull = zeros(nTraces, 1);

for k = 1:size(dataTraces,1)
    %% subtracts template waveform from data at different latencies
    testTrace = squeeze(dataTraces(k, centerChannel, :));
    [subtractedTraces offsets] = subtractWithShifts(testTrace, template, shiftStart, shiftEnd, shiftStep);
    
    %% fitting exponential to the artifact (traces after subtraction) 

    errors = zeros(1, length(subtractedTraces));
    for i = 1:length(subtractedTraces)
        errors(i) = norm(artifactMean(10:nSamples)-subtractedTraces{i}(10:nSamples));
    end
    
    optOffsetsFull(k) = offsets(find(errors == min(errors), 1));
    optSubtractedFull{k} = subtractedTraces{find(errors == min(errors), 1)};
    optErrorsFull(k) = min(errors);
end


%% checking similarity of subtracted traces to hypothetical artifact

distFromArtifactSubtracted = zeros(1, nTraces);
distFromArtifactOriginal = zeros(1, nTraces);
successes = zeros(1, nTraces);
actualArtifact = zeros(nSamples, 1); %average of failure traces, based on full algorithm


for i = 1:nTraces
    %distFromArtifactSubtracted(i) = norm(optSubtractedFull{i} - artifactMean);
    distFromArtifactSubtracted(i) = optErrorsFull(i);
    distFromArtifactOriginal(i)   = norm(squeeze(dataTraces(i, centerChannel,10:nSamples)) - artifactMean(10:nSamples));
    successes(i) = distFromArtifactSubtracted(i) < distFromArtifactOriginal(i) && optOffsetsFull(i) ~= shiftStart && optOffsetsFull(i) ~= shiftEnd;
    if successes(i)
        latencies(i) = optOffsetsFull(i) + templateMinPos;
    else
        actualArtifact = actualArtifact + squeeze(dataTraces(i, centerChannel, :));
    end
end
actualArtifact = actualArtifact/sum(~successes);

successes = cast(successes, 'logical'); %not sure why this is necessary, but get an error otherwise



%% plots successes and failures

if saveFigures == 1;
    if usePrevArtifact
        h = summaryPlotGenerator(dataTraces, template, centerChannel, actualArtifact,...
            actualArtifactExpOnly, latencies, latenciesInitial, patternNumber, movieNumber, successes,...
            successesExpFit, pathToEi, channelRadius, neuronID, dataPath, artifactMean, previousArtifact);
        saveas(h, ['p' num2str(patternNumber) 'm' num2str(movieNumber) '_' num2str(startArtFit) '-' num2str(endArtFit) '.pdf'])
    elseif artifactSubtraction
        h = summaryPlotGenerator(dataTraces, template, centerChannel, actualArtifact,...
            actualArtifactExpOnly, latencies, latenciesInitial, patternNumber, movieNumber, successes,...
            successesExpFit, pathToEi, channelRadius, neuronID, dataPath, artifactMean);
        saveas(h, ['p' num2str(patternNumber) 'm' num2str(movieNumber) '_' num2str(startLinFit) '-' num2str(endLinFit) '.pdf'])
    else
        h = summaryPlotGenerator(dataTraces, template, centerChannel, actualArtifact,...
            actualArtifactExpOnly, latencies, latenciesInitial, patternNumber, movieNumber, successes,...
            successesExpFit, pathToEi, channelRadius, neuronID, dataPath, artifactMean);
        saveas(h, ['p' num2str(patternNumber) 'm' num2str(movieNumber) '_' num2str(startExpFit) '-' num2str(endExpFit) '.pdf'])
    end
    close
end











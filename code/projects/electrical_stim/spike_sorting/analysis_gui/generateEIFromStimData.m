function ei = generateEIFromStimData(elecResp, movieIndex, varargin)
% to be used with addActiveNeuron



p = inputParser;

p.addRequired('elecResp', @isstruct)
p.addRequired('movieIndex', @isnumeric)


p.addParamValue('xLimits', [], @isnumeric)

p.parse(elecResp, movieIndex, varargin{:})

xLimits = p.Results.xLimits;


centerChannel = elecResp.cells.recElec;
currentMovie = elecResp.stimInfo.movieNos(movieIndex);

if isnumeric(elecResp.stimInfo.patternNo)
    patternNo = num2str(elecResp.stimInfo.patternNo);
else
    patternNo = elecResp.stimInfo.patternNo;
end

%loads data from the preprocessed data files
if exist([elecResp.names.data_path filesep 'p' patternNo filesep...
        'p' patternNo '_m' num2str(currentMovie)], 'file')
    dataTraces=NS_ReadPreprocessedData([elecResp.names.data_path filesep 'p' patternNo],...
        '', 0, patternNo, currentMovie, 99999);
elseif exist([elecResp.names.data_path filesep 'p' patternNo '_m' num2str(currentMovie)], 'file')
    dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, patternNo,...
        currentMovie, 99999);
else
    warnH = warndlg([elecResp.names.data_path filesep 'p' patternNo '_m' num2str(currentMovie) ' couldn''t be found']);
    ei = [];
    uiwait(warnH)
    return
end

nTraces = size(dataTraces, 1);
nElecs = size(dataTraces, 2);
nSamples = size(dataTraces, 3);

if length(xLimits) ~= 2
    xLimits = [1 nSamples];
end

actualArtifactFull = zeros(1, nElecs, nSamples); %average of failure traces, based on full algorithm
failures = 0;
for i = 1:nTraces
    if ~elecResp.analysis.latencies{movieIndex}(i)
        actualArtifactFull = actualArtifactFull + (dataTraces(i, :, :));
        failures = failures + 1;
    end
end
actualArtifactFull = actualArtifactFull/failures;



%subtracts estimated artifact
% if isempty(elecResp.analysis.estArtifact{movieIndex})
%     warnh = warndlg('No artifact estimate has been calculated yet.  Artifact estimate is required to generate a new EI.');
%     uiwait(warnh)
%     ei = [];
%     return
% else
subtractionVector = squeeze(actualArtifactFull(1, centerChannel, :));
% end

dataToPlot = zeros(nTraces, xLimits(2)-xLimits(1)+1);
for i = 1:elecResp.stimInfo.nPulses(movieIndex)
    dataToPlot(i, :) = squeeze(dataTraces(i, centerChannel, xLimits(1):xLimits(2))) - subtractionVector(xLimits(1):xLimits(2));
end

traceIndeces = chooseTracesGui(dataToPlot, 'returnSingleIndex', 'false');

traceMinIndsTemp = zeros(length(traceIndeces), 1);
traceMinInds = zeros(length(traceIndeces), 1);


for i = length(traceIndeces):-1:1
    traceMinIndsTemp(i) = find(dataToPlot(traceIndeces(i),:) == min(dataToPlot(traceIndeces(i),:)), 1);
    traceMinInds(i) = traceMinIndsTemp(i) + (xLimits(1)-1);
    if traceMinInds(i) <= 10 || traceMinInds(i) >= xLimits(2)-xLimits(1) + 1 - 15 %exclude traces that have minima too close to edges of trace
        traceMinInds(i) = [];
        traceMinIndsTemp(i) = [];
        traceIndeces(i) = [];
    end
end
    
isAligned = checkAlignment(dataToPlot(traceIndeces, :), traceMinIndsTemp);
if isAligned == 1
    %calculate the ei!
    %first subtract the estimated artifact
    subtractedTraces = zeros(length(traceIndeces), nElecs, nSamples);
    for i = 1:length(traceIndeces)
        subtractedTraces(i, :, :) = dataTraces(traceIndeces(i), :, :) - actualArtifactFull;
    end
    
    %create a matrix of aligned traces from which average (ei) will be calculated
    eiAllTraces = inf*ones(length(traceIndeces), nElecs, 81);
    for i = 1:length(traceMinInds)
        %determines relevant chunk of trace (places minimum at sample 19 of a 81-sample trace)
        firstTraceSample = traceMinInds(i)+1 - 19; %putting this sample in position 1 will put minimum in position 19
        lastTraceSample = firstTraceSample + 80;
        
        %determines where trace will be inserted into ei, so that "empty" regions of time are left
        %as infinities (flag for later averaging)
        if firstTraceSample < 1
            firstEiSample = 1 - firstTraceSample + 1;
            firstTraceSample = 1;
        else
            firstEiSample = 1;
        end
        
        if lastTraceSample > nSamples
            lastEiSample = 81 - (lastTraceSample - nSamples);
            lastTraceSample = nSamples;
        else
            lastEiSample = 81;
        end

        eiAllTraces(i, :, firstEiSample:lastEiSample) = subtractedTraces(i, :, firstTraceSample:lastTraceSample);
    end
    
    ei = zeros(nElecs, 81);
    for i = 1:81 %for each sample in ei, average traces that are not infinity (indicates that the region existed in the stim data)
        nInSum = 0;
        for j = 1:length(traceMinInds)
            if ~isinf(eiAllTraces(j, 1, i)) 
                nInSum = nInSum + 1;
                ei(:, i) = ei(:, i) + squeeze(eiAllTraces(j, :, i)');
            end
        end
        if nInSum > 0
            ei(:,i) = ei(:,i)/nInSum;
        else
            ei(:,i) = 0;
        end
    end
    

else
    ei = [];
    return
end
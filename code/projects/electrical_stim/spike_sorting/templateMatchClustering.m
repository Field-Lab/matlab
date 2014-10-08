function elecResp = templateMatchClustering(elecResp, movieNumber, varargin)

% params.shiftLimSpikeMin: limits of template offsets explored, in samples from beginning of data to
%   template minimum
% params.residLim: defines limits of error (2-norm of residual) calculation; start should be about
%   shiftStart + 20 (minimum 10), end should be about shiftEnd + 30 or larger
%
% latency values reflect position of spike minimum, where first sample of
% dataTrace is sample 1
%
%
% 2009-06: now able to use multiple electrodes to discriminate between spikes/artifact
% 2009-06: added "hand-picked" artifact model as a model type (user chooses one trace to serve as
% the initial model of the artifact)
%
% 2012-06-07: for some reason, traces with flags were being used to
% determine second artifact estimate (was commented out) - uncommented

%% argument parsing

p = inputParser;

p.addRequired('elecResp', @isstruct)
p.addRequired('movieNumber', @isnumeric)

p.addParamValue('saveFigures', false)
p.addParamValue('shiftLimSpikeMin', [10 40], @isnumeric)
p.addParamValue('shiftStep', 0.25, @isnumeric)
p.addParamValue('residLim', [11 40], @isnumeric);
p.addParamValue('modelType', 'linkage', @(x)any(strcmpi(x,{'linkage', 'prevArtifact', 'nextArtifact', 'currentArtifact', 'singExp', 'sumExp', 'ttx', 'handPicked', 'otherOccurrences'})))
p.addParamValue('progressBar', 'off', @(x)any(strcmpi(x,{'on', 'off', 'true', 'false'})))
p.addParamValue('analysisFlags',...
    elecResp.analysis.details.analysisFlags{movieNumber == elecResp.stimInfo.movieNos}, @(x)isnumeric(x))
p.addParamValue('noRefit', false, @islogical)
p.addParamValue('spikeMinExclusionReg', 10, @isnumeric) %number of samples from template offset boundary to exclude
% (prevent from being minimum error offset)  when initial minimum error
% occurs at one of the offset boundaries

p.parse(elecResp, movieNumber, varargin{:})

modelType = p.Results.modelType;
progressBar = p.Results.progressBar;
analysisFlags = p.Results.analysisFlags;
noRefit = p.Results.noRefit;
spikeMinExclusionReg = p.Results.spikeMinExclusionReg;

prevLatencies = [elecResp.analysis.latencies{movieNumber == elecResp.stimInfo.movieNos} elecResp.analysis.otherLatencies{movieNumber == elecResp.stimInfo.movieNos}];

if any(strcmpi(progressBar, {'on', 'true'}))
    dispProgressBar = 1;
else
    dispProgressBar = 0;
end

saveFigures = p.Results.saveFigures;
shiftLimSpikeMin = p.Results.shiftLimSpikeMin;
shiftStep = p.Results.shiftStep;
startFit = p.Results.residLim(1); 
endFit = p.Results.residLim(2);

% extracting info from elecResp
dataPath = elecResp.names.data_path;
patternNumber = elecResp.stimInfo.patternNo;
%pathToEi = elecResp.names.rrs_ei_path;
neuronID = [elecResp.cells.main elecResp.cells.active{movieNumber == elecResp.stimInfo.movieNos}];
centerChannel = elecResp.cells.recElec;
goodChannels = elecResp.cells.goodElecs;

% checks for previous elecResp format
if any(elecResp.cells.active{movieNumber == elecResp.stimInfo.movieNos} == elecResp.cells.main)
    errordlg('neuron appears in elecResp.cells.main and elecResp.cells.active')
end

centerChannelIndex = find(goodChannels == centerChannel);
nElecs = length(goodChannels);
channelRadius = 1;
makeSliderPlot = 0;

%% checks for user errors

movieIndex = find(movieNumber == elecResp.stimInfo.movieNos);

% checks to make sure that current artifact exists, otherwise defaults to using previous artifact
if (~isempty(elecResp.analysis.estArtifact{movieIndex}))...
        && strcmpi(modelType, 'currentArtifact') && ~isnan(max(max(elecResp.analysis.estArtifact{movieIndex})))
    estArtifact = elecResp.analysis.estArtifact{movieIndex};
elseif strcmpi(modelType, 'currentArtifact')
    modelType = 'prevArtifact';
end


% checks to make sure the previous artifact is valid, otherwise defaults to using linkage-based
% model
if (movieIndex > 1 && ~isempty(elecResp.analysis.estArtifact{movieIndex-1}))...
        && strcmpi(modelType, 'prevArtifact') && ~isnan(max(max(elecResp.analysis.estArtifact{movieIndex-1})))
    estArtifact = elecResp.analysis.estArtifact{movieIndex-1};
elseif strcmpi(modelType, 'prevArtifact')
    modelType = 'linkage';
end

% checks to make sure the next artifact is valid, otherwise defaults to using linkage-based
% model
if (movieIndex < length(elecResp.stimInfo.movieNos) && ~isempty(elecResp.analysis.estArtifact{movieIndex+1}))...
        && strcmpi(modelType, 'nextArtifact') && ~isnan(max(max(elecResp.analysis.estArtifact{movieIndex+1})))
    estArtifact = elecResp.analysis.estArtifact{movieIndex+1};
elseif strcmpi(modelType, 'nextArtifact')
    modelType = 'linkage';
end

    
% if any of the analysis flags is 1 (skip reanalysis), checks to make sure data has been previously
% analyzed
if any(any(analysisFlags==1)) && isempty(prevLatencies)
    warnH = warndlg(['Data must be previously analyzed (elecResp.analysis.latencies must not be empty)'...
        'in order to skip the reanalysis of some traces/templates']);
    uiwait(warnH)
    return
end

%% parameters for linkage analysis

linkageWindow = [10 40];
linkThresh = 30;
nBranches = 2;


%% loads data

%returns data traces: nTraces x nElectrodes x nSamples
dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, patternNumber, movieNumber, 99999);

nTraces  = size(dataTraces, 1);
nSamples = size(dataTraces, 3);
nTemplates = length(neuronID);

elecResp.stimInfo.nPulses(elecResp.stimInfo.movieNos == movieNumber) = nTraces;

if isempty(analysisFlags) || size(analysisFlags, 1) ~= elecResp.stimInfo.nPulses(elecResp.stimInfo.movieNos == movieNumber)...
        || size(analysisFlags, 2) ~= nTemplates
    analysisFlags = zeros(elecResp.stimInfo.nPulses(elecResp.stimInfo.movieNos == movieNumber), nTemplates);
    disp('analysis flags was reset to all zeros due to improper size')
end

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
% end
% clear eiFile

% calculates template shift limits based on minima (expands region for some neurons if they have
% different minima)

shiftStart = shiftLimSpikeMin(1) - max(templateMinPos);
shiftEnd = shiftLimSpikeMin(2) - min(templateMinPos);
nShifts = (shiftEnd - shiftStart)/shiftStep + 1;

%% calculated indices within shifts that correspond with latencies from previous analysis

prevOffsetInd = zeros(nTraces, nTemplates);
if any(any(analysisFlags == 1))
    offsets = shiftStart:shiftStep:shiftEnd;
    prevOffsetInd = zeros(nTraces, nTemplates);
    for i = 1:nTemplates
        for j = 1:nTraces
            if prevLatencies(j,i)
                prevOptOffset = prevLatencies(j,i) - templateMinPos(i);
                prevOffsetInd(j,i) = find(offsets == prevOptOffset);
            end
        end
    end
end


%%
if strcmpi(modelType, 'handPicked')
    
    dispTraces = [];
    for i = 1:nElecs
        dispTraces = [dispTraces squeeze(dataTraces(:, goodChannels(i), :))]; %#ok<AGROW>
    end
    
    selectedTraces = chooseTracesGui(dispTraces);
  
    estArtifact = reshape(mean(dataTraces(selectedTraces, goodChannels, :), 1), length(goodChannels), []);
    %figure
    %plot(estArtifact(1,:))
end

%%
if strcmpi(modelType, 'otherOccurrences')
    
    %find other elecResps with same pattern (other "occurances")
    oInd = strfind(elecResp.stimInfo.patternNo, 'o');
    
    if isempty(oInd)
        error('''other occurrences'' option is only available for elecResps labelled as occurrences')
    end
    
    files = dir([elecResp.names.data_path filesep 'elecResp*p' (elecResp.stimInfo.patternNo(1:oInd-1)) '*']);
    
    %collect potential traces from other occurances with finalized analysis of same amplitude
    %stimulus
    otherOTraces = [];
    
    thisAmp = elecResp.stimInfo.stimAmps(elecResp.stimInfo.movieNos == movieNumber);
    for ii = 1:length(files)
        %load other elecResp
        temp = load([elecResp.names.data_path filesep files(ii).name]);
        
        %find which movie has same amplitude (within 1%)
        matchInd = find(abs((temp.elecResp.stimInfo.stimAmps - thisAmp)/thisAmp) < 0.01);
        if length(matchInd) ~= 1
            error('One of the other occurrences doesn''t have a single matching stimulus amplitude')
        end
        
        if temp.elecResp.analysis.finalized(matchInd) && sum(temp.elecResp.analysis.latencies{matchInd} == 0) > 0 ...
                && ~strcmp(elecResp.stimInfo.patternNo, temp.elecResp.stimInfo.patternNo)
            
            dataTracesTemp = NS_ReadPreprocessedData(dataPath, '', 0, temp.elecResp.stimInfo.patternNo, temp.elecResp.stimInfo.movieNos(matchInd), 99999);
            
            artInd = find(temp.elecResp.analysis.latencies{matchInd} == 0);
            for jj = 1:length(artInd)
                otherOTraces = cat(1, otherOTraces, dataTracesTemp(artInd, :, :));
            end
            
        end
        
        clear dataTracesTemp temp
        
        %break loop if > 100 traces have been collected
        if size(otherOTraces, 1) > 100
            break
        end
    end
    

    if (movieIndex > 1 && ~isempty(elecResp.analysis.estArtifact{movieIndex-1})) && ~isnan(max(max(elecResp.analysis.estArtifact{movieIndex-1})))
        prevEstArtifact = elecResp.analysis.estArtifact{movieIndex-1};
    end
    
    
    otherOTracesToPlot = squeeze(otherOTraces(:, centerChannel, :));
    if exist('prevEstArtifact', 'var')
        selectedTraces = chooseTracesGui(otherOTracesToPlot, 'figTitle', 'remove outlier(s)', 'plotAlso', prevEstArtifact);
    else
        selectedTraces = chooseTracesGui(otherOTracesToPlot, 'figTitle', 'remove outlier(s)');    
    end
    otherOEstArtifact = squeeze(mean(otherOTraces(selectedTraces, :, :), 1));
    estArtifact = otherOEstArtifact(goodChannels, :);
end


%% TTX artifact subtraction

if strcmpi(modelType, 'ttx') %for now, assumes pattern number of data and ttx artifact correspond
    artifactTraces = NS_ReadPreprocessedData(elecResp.names.artifact_path, '', 0, patternNumber,...
        elecResp.stimInfo.artMovieNos(movieIndex), 99999);
    estArtifact = mean(artifactTraces, 1);
    estArtifact = reshape(estArtifact(1, goodChannels, :), length(goodChannels), []);
end


%% generates estimated artifact based on linkage analysis and subtracts it from data


if strcmpi(modelType, 'linkage')
    linkData = dataTraces(:, goodChannels, :);
    linkTemplate{1} = templates{1}; % for now, estArtWithLinkage can't handles multiple templates
    estArtifact = estArtWithLinkage(linkData, linkageWindow, linkThresh, nBranches, linkTemplate,...
        shiftStart, shiftEnd, shiftStep, centerChannelIndex);
end

%% loops through traces in data

% analysisFlags
% array of flag (same dimensions as latencies)
% = 0 means analyze normally
% = 1 means skip analysis of trace/template (leave as is)
% = 2 means force trace to not be failure for specific template
% = 3 means force trace to not be a failure for one of the possible templates (must = 3 for all
% templates)
% = 4 means force trace to be a failure for specific template


optOffsetsInit = zeros(nTraces, nTemplates);
optSubtracted = cell(1, nTraces);
optErrorsInit = zeros(1, nTraces);
allErrors = cell(1, nTraces);

if dispProgressBar
    progressBarH = figure('position', [500 600 500 50], 'Toolbar', 'none', 'Menubar', 'none',...
        'Color', 'white', 'Resize', 'off', 'Name', 'analyzing...');
    % 'WindowStyle', 'modal'
    progressAxes = axes('Parent', progressBarH, 'Units', 'pixels', 'position', [10 10 480 30],...
        'visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
    progressColors = cool(nTraces);
    hold on
end

for k = 1:nTraces
    
    if dispProgressBar
        try
            axes(progressAxes)
            current = plot([0 k/nTraces], [0.75 0.75], 'LineWidth', 10);
            set(findobj(current,'Type','line'), 'Color', progressColors(k,:))
        catch
        end
    else
        %disp(['computing for trace ' num2str(k)])
    end

    
    if ~all(analysisFlags(k,:) == 1)
        %% subtracts template waveform from data at different latencies
        testTrace = reshape(dataTraces(k, goodChannels, :), length(goodChannels), []);
        [subtractedTraces offsets] = subtractWithShifts(testTrace, templates, shiftStart, shiftEnd, shiftStep);

        
        %% finding error between model artifact and traces after subtraction
        if exist('estArtifact', 'var')
            errors = calcErrors(subtractedTraces, modelType, startFit, endFit, estArtifact,...
                analysisFlags(k,:), prevOffsetInd(k,:));
        else
            errors = calcErrors(subtractedTraces, modelType, startFit, endFit, [],...
                analysisFlags(k,:), prevOffsetInd(k,:));
        end
        
        allErrors{k} = errors;
        
        optErrorsInit(k) = min(reshape(errors,[],1));
        minInd = find(errors == optErrorsInit(k),1);
        
        
        if makeSliderPlot
            % makes slider plot (only works for one electrode, one neuron cases)
            sliderFig = figure;
            slider = make_loop_slider_list(1,1,length(subtractedTraces));

            if size(testTrace, 1) > 1
                plotTestTrace = squeeze(testTrace(1,:));
                plotSubtractedTraces = cell(length(subtractedTraces), 1);
                for i = 1:length(subtractedTraces)
                    plotSubtractedTraces{i} = squeeze(subtractedTraces{i}(1,:));
                end
                plotEstArtifact = squeeze(estArtifact(1,:));
            else
                plotTestTrace = testTrace;
                plotSubtractedTraces = subtractedTraces;
                plotEstArtifact = estArtifact;
            end
            
            while ishandle(sliderFig)
                i = round(get(slider,'Value'));
                cla
                hold on

                plot(plotTestTrace, 'k')
                plot(plotSubtractedTraces{i}, 'm')
                plot(plotTestTrace - plotSubtractedTraces{i}, 'b')
                plot(plotEstArtifact, 'g')
                plot(plotTestTrace - plotEstArtifact, 'r')
                legend('raw data', 'subtracted', 'template', 'estimated artifact', '')
                title(['offset number ', num2str(i-1), ', error = ', num2str(errors(i)), 10, ...
                    'lowest error occurs at offset number ', num2str(minInd-1)])
                
                plot([startFit startFit], [min(plotTestTrace) max(plotTestTrace-plotSubtractedTraces{i})], 'k--')
                plot([endFit endFit], [min(plotTestTrace) max(plotTestTrace-plotSubtractedTraces{i})], 'k--')
                hold off
                uiwait;
            end
            keyboard % this is your opportunity to set makeSliderPlot to 0!!!
        end


        % a horrible kluge to make ind2sub actually WORK
        if nTemplates < 2
            minSub(1) = find(errors == optErrorsInit(k), 1);
        else
            if nTemplates < 3
                [minSub(1) minSub(2)] = ind2sub(size(errors), find(errors == optErrorsInit(k), 1));
            else
                if nTemplates < 4
                    [minSub(1) minSub(2) minSub(3)] = ind2sub(size(errors), find(errors == optErrorsInit(k), 1));
                else
                    [minSub(1) minSub(2) minSub(3) minSub(4)] = ind2sub(size(errors), find(errors == optErrorsInit(k), 1));
                end
            end
        end

        % removing error minima that lie at smallest and largest offset values
        noMoreExtremeMins = 0;

        while ~noMoreExtremeMins
            if any(minSub==2) %if lowest offset gives lowest error, set it and the offsets within the next spikeMinExclusionReg samples of that template to infinity
                errors(minInd) = inf;
                badDim = find(minSub==2, 1);
                for i = 1:(spikeMinExclusionReg/shiftStep)
                    afterMinInd = minInd + i*size(errors,1)^(badDim-1);
                    errors(afterMinInd) = inf;
                end
            elseif any(minSub==(nShifts+1)) %if highest offset gives lowest error, set it and the offsets within the previous spikeMinExclusionReg samples of that template to infinity
                errors(minInd) = inf;
                badDim = find(minSub == nShifts+1, 1);
                for i = 1:(spikeMinExclusionReg/shiftStep)
                    beforeMinInd = minInd - i*size(errors,1)^(badDim-1);
                    errors(beforeMinInd) = inf;
                end
            else
                noMoreExtremeMins = 1;
            end

            optErrorsInit(k) = min(reshape(errors,[],1));
            minInd = find(errors == optErrorsInit(k),1);

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
            optOffsetsInit(k,i) = offsets(minSub(i));
        end
        optSubtracted{k} = subtractedTraces{minInd};
    else  %use latency values from previous analysis
        for i = 1:nTemplates
            optOffsetsInit(k,i) = prevLatencies(k,i) - templateMinPos(i); 
        end
    end
end


%% chooses error that is lower from optErrorsInit and errorsNoSub

%%%%%% use only non-flagged, artifact-alone traces to generate updated
%%%%%% artifact estimate; if none are available default to previous version
%%%%%% (include success traces with spike templates subtracted)

successesInit = zeros(nTraces, nTemplates);
actualArtifactInit = zeros(nElecs, nSamples); %keep this separate from artifact estimate for refitting just to be consistent with how final artifact estimate is calculated
artForRefit = zeros(nElecs, nSamples);
latenciesInit = zeros(nTraces, nTemplates);

failures = 0;
nArtTraces = 0;
for i = 1:nTraces
    if all(isnan(optOffsetsInit(i,:)))
        failures = failures + 1;
        actualArtifactInit = actualArtifactInit + reshape(dataTraces(i, goodChannels, :), length(goodChannels), []);
        if all(analysisFlags(i,:) == 0)
            nArtTraces = nArtTraces+1;
            artForRefit = artForRefit + reshape(dataTraces(i, goodChannels, :), length(goodChannels), []);
        end
    else
        for j = 1:nTemplates            
            if ~isnan(optOffsetsInit(i,j))
                latenciesInit(i,j) = optOffsetsInit(i,j) + templateMinPos(j);
                successesInit(i,j) = 1;
            end
        end
    end
end

actualArtifactInit = actualArtifactInit/failures;
artForRefit = artForRefit/nArtTraces;
successesInit = cast(successesInit, 'logical'); %not sure why this is necessary, but get an error otherwise


%if ~any(any(isnan(actualArtifactFull))) && max(max(abs(actualArtifactFull))) %only use actual artifact if it could be calculated (>=1 failure trace)

%% creates hypothetical artifact from best fits to estimated artifact (from subtracted and non-subtracted)
%%% previously (up to June 2012) this artifact estimate was used for
%%% refitting in all cases, but now it's only used if there are no traces
%%% that are both 1) failures for all analyzed neurons and 2) unflagged


if ~noRefit
  
    if any(any(isnan(artForRefit))) || max(max(abs(artForRefit)))==0 %no unflagged failure traces were available to generate artifact estimate for refitting
        %optErrorsInitInclude = optErrorsInit(sum(analysisFlags, 2) == 0); %don't calculate mean using errors from traces that are not analyzed normally
        %meanOptError = mean(optErrorsInitInclude);
        artifactSum = zeros(nElecs, nSamples);
        
        %collapse flags across cells being analyzed
        anyFlags = any(analysisFlags,2);
                
        nArtifacts = 0;
        for i = 1:nTraces
            if ~anyFlags(i) || all(anyFlags)%added this criterion June 2012, only used if there are unflagged traces
                artifactSum = artifactSum + optSubtracted{i};
                nArtifacts = nArtifacts + 1;
            end
        end
        
        artForRefit = artifactSum / nArtifacts;
        %artifactMean = artifactSum / nArtifacts;
        
        clear artifactSum nArtifacts %optErrorsInit optErrorsInitInclude meanInitError
    end
    %% redoing error minimization of template subtraction at different latencies based on hyp. artifact
    
    optOffsetsFull = zeros(nTraces, nTemplates);
    
    optSubtractedFull = cell(nTraces, 1);
    optErrorsFull = zeros(nTraces, 1);
    
    for k = 1:nTraces
        
        if dispProgressBar
            try
                axes(progressAxes)
                current = plot([0 k/nTraces], [0.25 0.25], 'r', 'LineWidth', 10);
                set(findobj(current,'Type','line'), 'Color', progressColors(k,:))
            catch
            end
        end
        
        if ~all(analysisFlags(k,:) == 1)
            
            testTrace = reshape(dataTraces(k, goodChannels, :), length(goodChannels), []);
            [subtractedTraces offsets] = subtractWithShifts(testTrace, templates, shiftStart, shiftEnd, shiftStep);
            
            errors = calcErrors(subtractedTraces, 'hypArtifact', startFit, endFit, artForRefit,...
                analysisFlags(k,:), prevOffsetInd(k,:));
            
            optErrorsFull(k) = min(reshape(errors,[],1));
            minInd = find(errors == optErrorsFull(k),1);
            
            % a horrible kluge to make ind2sub actually WORK
            if nTemplates < 2
                minSub(1) = find(errors == optErrorsFull(k), 1);
            else
                if nTemplates < 3
                    [minSub(1) minSub(2)] = ind2sub(size(errors), find(errors == optErrorsFull(k), 1));
                else
                    if nTemplates < 4
                        [minSub(1) minSub(2) minSub(3)] = ind2sub(size(errors), find(errors == optErrorsFull(k), 1));
                    else
                        [minSub(1) minSub(2) minSub(3) minSub(4)] = ind2sub(size(errors), find(errors == optErrorsFull(k), 1));
                    end
                end
            end
            
            % removing error minima that lie at smallest and largest offset values
            noMoreExtremeMins = 0;
            while ~noMoreExtremeMins
                if any(minSub==2) %if lowest offset gives lowest error, set it and the offsets within the next spikeMinExclusionReg samples of that template to infinity
                    errors(minInd) = inf;
                    badDim = find(minSub==2, 1);
                    for i = 1:(spikeMinExclusionReg/shiftStep)
                        afterMinInd = minInd + i*size(errors,1)^(badDim-1);
                        errors(afterMinInd) = inf;
                    end
                elseif any(minSub==(nShifts+1)) %if highest offset gives lowest error, set it and the offsets within the previous spikeMinExclusionReg samples of that template to infinity
                    errors(minInd) = inf;
                    badDim = find(minSub == nShifts+1, 1);
                    for i = 1:(spikeMinExclusionReg/shiftStep)
                        beforeMinInd = minInd - i*size(errors,1)^(badDim-1);
                        errors(beforeMinInd) = inf;
                    end
                else
                    noMoreExtremeMins = 1;
                end
                
                optErrorsFull(k) = min(reshape(errors,[],1));
                minInd = find(errors == optErrorsFull(k),1);
                
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
                optOffsetsFull(k,i) = offsets(minSub(i));
            end
            optSubtractedFull{k} = subtractedTraces{minInd};
        else  %use latency values from previous analysis
            for i = 1:nTemplates
                optOffsetsFull(k,i) = prevLatencies(k,i) - templateMinPos(i);
            end
        end
    end

end

%% calculates artifact based on final classification

if ~noRefit
    successesFull = zeros(nTraces, nTemplates);
    actualArtifactFull = zeros(nElecs, nSamples); %average of failure traces, based on full algorithm
    latenciesFull = zeros(nTraces, nTemplates);
    
    failures = 0;
    for i = 1:nTraces
        if all(isnan(optOffsetsFull(i,:)))
            actualArtifactFull = actualArtifactFull + reshape(dataTraces(i, goodChannels, :), length(goodChannels), []);
            failures = failures + 1;
        else
            for j = 1:nTemplates
                if ~isnan(optOffsetsFull(i,j))
                    latenciesFull(i,j) = optOffsetsFull(i,j) + templateMinPos(j);
                    successesFull(i,j) = 1;
                end
            end
        end
    end
    
    actualArtifactFull = actualArtifactFull/failures;
    successesFull = cast(successesFull, 'logical'); %not sure why this is necessary, but get an error otherwise
else %just use initial classification
    successesFull = successesInit;
    actualArtifactFull = actualArtifactInit; %artifact estimate based on first round of classification
    latenciesFull = latenciesInit;
    artForRefit = estArtifact; %initial artifact estimate (pre-classification)
end
successRate = sum(successesFull(:,1))/nTraces;


%% updating elecResp

elecResp.analysis.type{movieIndex} = modelType;
elecResp.analysis.latencies{movieIndex} = latenciesFull(:,1);
if nTemplates > 1
    elecResp.analysis.otherLatencies{movieIndex} = latenciesFull(:,2:end);
end
elecResp.analysis.successRates(movieIndex) = successRate;
elecResp.analysis.erfCurrent = 0;
elecResp.analysis.finalized(movieIndex) = 0;

elecResp.analysis.details.residCalcWindow{movieIndex} = [startFit endFit];
elecResp.analysis.details.tempOffsetWindow{movieIndex} = [shiftStart shiftEnd];
elecResp.analysis.details.tempOffsetStep{movieIndex} = shiftStep;

elecResp.analysis.details.analysisFlags{movieIndex} = analysisFlags;

if ~any(any(isnan(actualArtifactFull))) && max(max(abs(actualArtifactFull))) %only use actual artifact if it could be calculated (>=1 failure trace)
    elecResp.analysis.estArtifact{movieIndex} = actualArtifactFull;
else
    elecResp.analysis.estArtifact{movieIndex} = artForRefit;
end


% analysis type-specific

if strcmpi(modelType, 'linkage')
    elecResp.analysis.details.linkCalcWindow{movieIndex} = linkageWindow;
    elecResp.analysis.details.linkThresh{movieIndex} = linkThresh;
    elecResp.analysis.details.nBranches{movieIndex} = nBranches;
else
    elecResp.analysis.details.linkCalcWindow{movieIndex} = [];
    elecResp.analysis.details.linkThresh{movieIndex} = [];
    elecResp.analysis.details.nBranches{movieIndex} = [];
end

if strcmpi(modelType, 'handPicked')
    elecResp.analysis.details.handPickedIndeces{movieIndex} = selectedTraces;
else
    elecResp.analysis.details.handPickedIndeces{movieIndex} = [];
end

if strcmpi(modelType, 'otherOccurrences')
    elecResp.analysis.details.otherOArtifact{movieIndex} = otherOEstArtifact;
else
    elecResp.analysis.details.otherOArtifact{movieIndex} = [];
end


%% plots successes and failures

% packaging variables into structures
if saveFigures
    artifact.actual     = actualArtifactFull;
    artifact.actualInit = actualArtifactInit;
    artifact.mean       = artForRefit;

    results.latenciesInit = latenciesInit;
    results.latenciesFull = latenciesFull;
    results.successesInit = successesInit;
    results.successesFull = successesFull;
    
    if any(strcmpi(modelType, {'prevArtifact', 'nextArtifact', 'linkage', 'ttx', 'hand-picked', 'currentArtifact'}))
        artifact.estimate = estArtifact;
    end
    
    h = summaryPlotGenerator(elecResp, dataTraces, channelRadius, movieNumber, artifact, templates, results);

    print(h, '-append', '-dpsc2', ['p', num2str(patternNumber), '_', num2str(startFit) '-' num2str(endFit), '.eps'])
    close
end

if dispProgressBar
    try
        close(progressBarH)
    catch
    end
end

end






function axesH = frequencyResponseSpikeWaveforms(elecRespPath, freq, movieInd, nSequence, elec, varargin)

% generates a raster plot of responses to pulse train at particular
% frequency
%
% ****note that 'trials' aren't really meaningful for frequencies that result in
% no gap in time between pulse sequences (e.g. 20 pulses at 20 Hz, repeated
% each second)
%
%
% a good portion of this code is dedicated to dealing with spikes from
% other cells that may contaminate the extracted waveforms (either by
% excluding these traces or subtracting the other spike waveform from the
% traces)
%


p = inputParser;

p.addRequired('elecRespPath', @ischar)
p.addRequired('freq', @isnumeric) %frequency of stimulus in Hz
p.addRequired('movieInd', @isnumeric)
p.addRequired('nSequence', @isnumeric) %number of pulses in a single repetition
p.addRequired('elec', @isnumeric)

p.addParamValue('otherSpikePath', '', @ischar) %path to elecResp containing another cell that responds, so that spikes occuring in these trials are dealt with
p.addParamValue('comparisonResponse', [], @isstruct) %contains .path, .movieInd, .freq, .otherSpikePath(optional) -- to use as reference spikes waveform
p.addParamValue('excludeOrSubtract', 'exclude', @(x)any(strcmpi(x,{'exclude', 'subtract', 'ignore'}))) %how to deal with spikes from cell in 'otherSpikePath'
p.addParamValue('traceLengthPlot', 30, @isnumeric) %how much of the each trace to plot
p.addParamValue('plotVsLatency', false, @islogical) %plot vs. spikes latency rather than vs. pulse number
p.addParamValue('yPlotOffset', 300, @isnumeric) %distance between different spikes
p.addParamValue('nthPulseToCompare', 0, @isnumeric) %mean uncontaminated spike waveforms corresponding to these pulses are plotted for comparison
                                                    %(0 means don't plot any, -1 means plot waveform from comparison frequency)

p.addParamValue('subtractPulseByPulseArtifact', false, @islogical)

p.parse(elecRespPath, freq, movieInd, nSequence, elec, varargin{:})


otherSpikePath = p.Results.otherSpikePath;
comparisonResponse = p.Results.comparisonResponse;

otherSpikeTreatment = p.Results.excludeOrSubtract;
plotVsLat = p.Results.plotVsLatency;

traceLengthPlot = p.Results.traceLengthPlot;
yPlotOffset = p.Results.yPlotOffset;

subtractPulseByPulseArtifact = p.Results.subtractPulseByPulseArtifact;
nthPulseToCompare = p.Results.nthPulseToCompare;


%%




%% calculated artifact-subtracted waveforms for test traces
    
[latencies successes] = extractFrequencyAnalysis(elecRespPath, movieInd, nSequence);
nReps = size(latencies,2);

% which trials have spikes from other cell (for exclusion or subtraction)
if ~isempty(otherSpikePath)
    [otherSpikeLat otherSpikeBin] = extractFrequencyAnalysis(otherSpikePath, movieInd, nSequence);
else
    otherSpikeBin = zeros(nSequence, nReps);
end

% load up the elecResp file and traces
load(elecRespPath)
dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
    elecResp.stimInfo.movieNos(movieInd), 99999);
dataTraces = squeeze(dataTraces(:,elec,:));

nTraces = size(dataTraces, 1); %traces are ordered according the actual time the pattern is played (all pulses in sequence for rep 1, then all pulses for rep 2, etc)
nSamples = size(dataTraces, 2);


% estimate artifact on chosen electrode and subtract from data (based on
% traces without target or nontarget cells' spikes)
failureBin = elecResp.analysis.latencies{movieInd} == 0 & ~reshape(otherSpikeBin,[],1);

if sum(failureBin)>0 %there are some traces without spikes, so just take mean of these
    failures = squeeze(dataTraces(failureBin, :));
    artifactEst = mean(failures,1);
elseif elec == elecResp.cells.recElec %use artifact estimate from elecResp because no spike-free traces are available
    artifactEst = elecResp.analysis.estArtifact{movieInd};
    warning(['no traces without spikes; using artifact estimate from elecResp file',...
        10, 'These may be contaminated by spikes from another cell!!']); %#ok<WNTAG>
else
    error(['no traces without spikes and no artifact estimate '...
    '(recording electrode used in elecResp is different from electrode chosen for this analysis)--aborting!']);
end
    


% subtract estimated artifact from all traces
subtractedData = zeros(nTraces,nSamples);

for ii = 1:nTraces
    subtractedData(ii,:) = dataTraces(ii,:) - artifactEst;
end



% hold on to subtracted artifact-only traces to look for changes in artifact
% waveform with pulse #, and subtract pulse-by-pulse artifact from data
% traces if specified in subtractPulseByPulseArtifact

if subtractPulseByPulseArtifact
    subtractedDataPulseByPulse = zeros(size(subtractedData));
end

artifactOnly = cell(1,nSequence);
for ii = 1:nSequence %pulses in sequence
    artifactOnly{ii} = [];
    for jj = 1:nReps; %sequence repetitions
        if failureBin((jj-1)*nSequence+ii) %if this traces doesn't have spike from target or nontarget cell
            artifactOnly{ii} = [artifactOnly{ii}; squeeze(subtractedData((jj-1)*nSequence+ii,:))];
        end
    end
    
    
    if subtractPulseByPulseArtifact
        if ~isempty(artifactOnly{ii})
            %figure
            %hold on
            %plot(mean(artifactOnly{ii},1), 'r-')
            for jj = 1:nReps
                subtractedDataPulseByPulse((jj-1)*nSequence+ii, :) = subtractedData((jj-1)*nSequence+ii, :) - mean(artifactOnly{ii},1);
                if successes(ii,jj)
             %       plot(subtractedData((jj-1)*nSequence+ii, :), 'k-')
             %       plot(subtractedDataPulseByPulse((jj-1)*nSequence+ii, :))
                end
            end
        end
    end
    
end



%% subtract out 'other' spikes if specified by otherSpikeTreatment

%figure
%hold on

%excludeSpikes = strcmpi(otherSpikeTreatment, 'exclude');
%includeSpikesNoSubtract = strcmp(otherSpikeTreatment, 'ignore');

if strcmpi(otherSpikeTreatment, 'subtract') && sum(sum(otherSpikeBin))>0
%if ~excludeSpikes && ~includeSpikesNoSubtract && sum(sum(otherSpikeBin))>0
    excludeElecResp = load(otherSpikePath);
    exRecElec = excludeElecResp.elecResp.cells.recElec;
    minPos = find(excludeElecResp.elecResp.cells.mainEI(exRecElec,:) ==...
        min(excludeElecResp.elecResp.cells.mainEI(exRecElec,:)));
    
    %get ei waveform for other neuron on specified electrode
    excludeWaveform{1} = excludeElecResp.elecResp.cells.mainEI(elec,:);

    clear excludeElecResp
    
    for ii = 1:nTraces
        offset = otherSpikeLat(ii) - minPos;
        if elecResp.analysis.latencies{movieInd}(ii) && otherSpikeBin(ii)
            %plot((1:50)+5*ii, subtractedData(ii,1:50), 'r')
            subtractedData(ii,:) = subtractWithShifts(subtractedData(ii,:), excludeWaveform, offset);
            if subtractPulseByPulseArtifact
                if any(subtractedDataPulseByPulse(ii,:)>0) %traces that don't have a valid artifact are set to zeros, so leave them as zeros
                    subtractedDataPulseByPulse(ii,:) = subtractWithShifts(subtractedDataPulseByPulse(ii,:), excludeWaveform, offset);
                end
            end
        end
        
%         if elecResp.analysis.latencies{movieInd}(ii) && otherSpikeBin(ii)
%             plot((1:50)+5*ii, subtractedData(ii,1:50), 'b')
%         elseif elecResp.analysis.latencies{movieInd}(ii)
%             plot((1:50)+5*ii, subtractedData(ii,1:50), 'k')
%         end
    end
end
%hold off
%% calculate artifact-subtracted waveform from comparison frequency

if ~isempty(comparisonResponse)

    temp = load(comparisonResponse.path);
    
    elecRespCompare = temp.elecResp;
    
    if isfield(comparisonResponse, 'otherSpikePath') && ~isempty(comparisonResponse.otherSpikePath)
        [x otherSpikeBinCompare] = extractFrequencyAnalysis(comparisonResponse.otherSpikePath, comparisonResponse.movieInd);
    else
        otherSpikeBinCompare = zeros(1, elecRespCompare.stimInfo.nPulses(comparisonResponse.movieInd));
    end
    
    dataTracesCompare=NS_ReadPreprocessedData(elecRespCompare.names.data_path, '', 0, elecRespCompare.stimInfo.patternNo,...
        elecRespCompare.stimInfo.movieNos(comparisonResponse.movieInd), 99999);

    %estimate artifact on chosen electrode and subtract from data
    failureBinCompare = elecRespCompare.analysis.latencies{comparisonResponse.movieInd} == 0 & ~otherSpikeBinCompare';
    failuresCompare = squeeze(dataTracesCompare(failureBinCompare, elec, :));    
    %failuresCompare = squeeze(dataTracesCompare(elecRespCompare.analysis.latencies{comparisonResponse.movieInd} == 0, elec, :));
    artifactEstCompare = mean(failuresCompare,1);
    subtractedDataCompare = zeros(0,size(dataTracesCompare,3));
    
    %figure
    %hold on
    for ii = 1:size(dataTracesCompare,1)
        if ~otherSpikeBinCompare(ii) && elecRespCompare.analysis.latencies{comparisonResponse.movieInd}(ii)
            subtractedDataCompare = [subtractedDataCompare; squeeze(dataTracesCompare(ii,elec,:))' - artifactEstCompare]; %#ok<AGROW>
            %plot(subtractedDataCompare(end,:))
        end
    end
    %hold off
end

%% plot (unexcluded) successes with position as a function of pulse # in
% sequence

%remove first 10 samples of each trace (contaminated by artifact)
subtractedData = subtractedData(:,11:end);
if subtractPulseByPulseArtifact
    subtractedDataPulseByPulse = subtractedDataPulseByPulse(:,11:end);
end
if exist('subtractedDataCompare', 'var')
    subtractedDataCompare = subtractedDataCompare(:,11:end);
end
for kk = 1:nSequence
    artifactOnly{kk} = artifactOnly{kk}(:,11:end);
end

if ~plotVsLat
    figure('position', [300 200 400 200])
else
    figure('position', [300 200 400 100])
end
axesH = axes;

hold on

% determine mean (uncontaminated) spike waveform from nth pulse in
% sequence

if any(nthPulseToCompare>0)
    compareWaveSum = zeros(1, size(subtractedData,2));
    nInSum = 0;
    for kk = nthPulseToCompare
        for jj = 1:nReps
            if successes(kk,jj) && ~otherSpikeBin(kk,jj)
                if ~subtractPulseByPulseArtifact
                    compareWaveSum = compareWaveSum + subtractedData((jj-1)*nSequence + kk, :);
                    nInSum = nInSum+1;
                else
                    compareWaveSum = compareWaveSum + subtractedDataPulseByPulse((jj-1)*nSequence + kk, :);
                    nInSum = nInSum+1;
                end
            end
        end
    end
    if nInSum == 0
        keyboard %there aren't any uncontaminated spikes from any of the pulses you want to compare with
    end
    compareWaveMean = compareWaveSum/nInSum;
elseif nthPulseToCompare == -1 %use comparison frequence spike waveform
    compareWaveMean = mean(subtractedDataCompare,1);
end

wavesPlotted = cell(nSequence);
minOffset = [];
for kk = 1:nSequence;
    
    wavesPlotted{kk} = [];
    for jj = 1:nReps
        if successes(kk,jj) %target cell spike!
            if ~(otherSpikeBin(kk,jj) && strcmpi(otherSpikeTreatment, 'exclude'))
                if ~subtractPulseByPulseArtifact
                    wavesPlotted{kk} = [wavesPlotted{kk}; subtractedData((jj-1)*nSequence + kk, :)]; %store for mean, std plotting
                else
                    wavesPlotted{kk} = [wavesPlotted{kk}; subtractedDataPulseByPulse((jj-1)*nSequence + kk, :)]; %store for mean, std plotting
                end
            end
            
            if ~otherSpikeBin(kk,jj) %uncomtaminated spike
                if ~plotVsLat && ~subtractPulseByPulseArtifact
                    plot((1:traceLengthPlot) + traceLengthPlot*kk, subtractedData((jj-1)*nSequence + kk, 1:traceLengthPlot), 'k')
                elseif subtractPulseByPulseArtifact
                    plot((1:traceLengthPlot) + traceLengthPlot*kk, subtractedDataPulseByPulse((jj-1)*nSequence + kk, 1:traceLengthPlot), 'k')
                else %plot vs. rounded latency instead of vs. pulse number
                    plot((1:traceLengthPlot)+traceLengthPlot*floor(latencies(kk,jj)), subtractedData((jj-1)*nSequence + kk, 1:traceLengthPlot), 'k')
                end
                    
                if isempty(minOffset) %hold on to for later plot specs
                    minOffset = find(subtractedData((jj-1)*nSequence + kk, :) == min(subtractedData((jj-1)*nSequence + kk, :)));
                end
                
            elseif otherSpikeBin(kk,jj) && ~strcmpi(otherSpikeTreatment, 'exclude') %subtracted or contaminated spike
                if ~plotVsLat && ~subtractPulseByPulseArtifact
                    plot((1:traceLengthPlot) + traceLengthPlot*kk, subtractedData((jj-1)*nSequence + kk, 1:traceLengthPlot), 'r')
                elseif subtractPulseByPulseArtifact
                    plot((1:traceLengthPlot) + traceLengthPlot*kk, subtractedDataPulseByPulse((jj-1)*nSequence + kk, 1:traceLengthPlot), 'r')
                else %plot vs. rounded latency instead of vs. pulse number
                    plot((1:traceLengthPlot)+traceLengthPlot*floor(latencies(kk,jj)), subtractedData((jj-1)*nSequence + kk, 1:traceLengthPlot), 'r')
                end
            end
        end
    end
    %plot mean of spikes from nth pulse in sequence for comparison
    if any(nthPulseToCompare)
        plot((1:traceLengthPlot) + traceLengthPlot*kk, compareWaveMean(1:traceLengthPlot), 'color', [0 0 1])
    end
end
% if all traces have subtracted spike waveforms, just use one of them to
% estimate spike min position
if isempty(minOffset)
    tempInd = find(successes,1);
    minOffset = find(subtractedData(tempInd,:) == min(subtractedData(tempInd,:)));
end


%determine and plot means and SDs of spike waveforms
%waveMeans = zeros(nSequence, size(wavesPlotted{1}, 2));
%waveSDs = zeros(nSequence, size(wavesPlotted{1}, 2));
if ~plotVsLat
    for kk = 1:nSequence
        if ~isempty(wavesPlotted{kk})
            %waveMeans(kk,:) = mean(wavesPlotted{kk},1);
            %waveSDs(kk,:) = std(wavesPlotted{kk},0,1);
            waveMean = mean(wavesPlotted{kk}(:,1:traceLengthPlot),1);
            waveSD = std(wavesPlotted{kk}(:,1:traceLengthPlot),0,1);
            
            xToPlot = (1:traceLengthPlot)+traceLengthPlot*kk;
            patch([xToPlot xToPlot(end:-1:1)],...
                [waveMean+waveSD waveMean(end:-1:1)-waveSD(end:-1:1)] + yPlotOffset, [0.75 0.75 0.75], 'edgeColor', 'none')
            plot(xToPlot, waveMean+yPlotOffset, 'k')
        end
        %plot mean of spikes from nth pulse in sequence for comparison
        if any(nthPulseToCompare)
            plot(xToPlot, compareWaveMean(1:traceLengthPlot)+yPlotOffset, 'color', [0 0 1])
        end
    end
end


% plot comparison spike waveforms, means and SDs
if exist('subtractedDataCompare', 'var')
    if ~plotVsLat
        xOffset = (nSequence+2)*traceLengthPlot;
    else
        xOffset = (max(max(latencies))+2)*traceLengthPlot;
    end
    
    for jj = 1:size(subtractedDataCompare, 1)
        plot((1:traceLengthPlot) + xOffset, subtractedDataCompare(jj,1:traceLengthPlot), 'k-')
    end
    
    if ~plotVsLat
        compMean = mean(subtractedDataCompare(:,1:traceLengthPlot),1);
        compSD = std(subtractedDataCompare(:,1:traceLengthPlot),0,1);

        xToPlot = (1:traceLengthPlot) + xOffset;
        patch([xToPlot xToPlot(end:-1:1)], [compMean+compSD compMean(end:-1:1)-compSD(end:-1:1)] + yPlotOffset, [0.75 0.75 0.75], 'edgeColor', 'none')
        plot(xToPlot, compMean+yPlotOffset, 'k')
    end
end

%plot artifact-only traces, and means +/- SD of artifact-only traces
if ~plotVsLat
    for kk = 1:nSequence
        if ~isempty(artifactOnly{kk})
            
            for jj = 1:size(artifactOnly{kk},1)
                plot((1:traceLengthPlot) + traceLengthPlot*kk, artifactOnly{kk}(jj,1:traceLengthPlot)-2*yPlotOffset, 'k-')
            end
            artMean = mean(artifactOnly{kk}(:,1:traceLengthPlot),1);
            artSD = std(artifactOnly{kk}(:,1:traceLengthPlot),0,1);

            xToPlot = (1:traceLengthPlot) + traceLengthPlot*kk;
            patch([xToPlot xToPlot(end:-1:1)], [artMean+artSD artMean(end:-1:1)-artSD(end:-1:1)] - yPlotOffset, [0.75 0.75 0.75], 'edgeColor', 'none')
            plot(xToPlot, artMean - yPlotOffset, 'k')
            
            text((nSequence+2)*traceLengthPlot, -2*yPlotOffset, ['residual', 10 'artifacts'])
        else
            plot(0, -2*yPlotOffset)
        end
    end
end


title(['responses for ' num2str(elecResp.stimInfo.stimAmps(movieInd)) 'uA, ' num2str(freq), ' Hz (' num2str(comparisonResponse.freq) ' Hz)'])

if ~plotVsLat
    xlabel('pulse number')
    set(gca, 'xtick', [31 nSequence*traceLengthPlot+1] + minOffset, 'xticklabel', [1 nSequence])
else
    xlabel('latency')
    set(gca, 'xtick', [], 'xticklabel', [], 'xlim', [(min(min(latencies(latencies~=0)))-1)*traceLengthPlot (max(max(latencies))+4)*traceLengthPlot])
end

hold off






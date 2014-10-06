%% example plots for frequency-dependent responses

clear all

% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n947_p58';
% elec = 60;
% elecRespPathBaseOtherSpikes = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n2_p58';
% 
% compareFreq = 10;
% comparisonResponse.movieInd = 10;
% nSequence = 20;

% 
% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n2_p58';
% elec = 64;
% elecRespPathBaseOtherSpikes = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n947_p58';
% 
% compareFreq = 10;
% comparisonResponse.movieInd = 15;


% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006/elecResp_n646_p36';
% elec = 26;
% 
% elecRespPathBaseOtherSpikes = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006/elecResp_n887_p36';
% compareFreq = 20;
% comparisonResponse.movieInd = 5;
% nSequence = 20;

% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028/elecResp_n875_p61';
% elec = 59;
% 
% elecRespPathBaseOtherSpikes = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028/elecResp_n2dummy_p61';
% compareFreq = 10;
% comparisonResponse.movieInd = 7;
% 
% nSequence = 20;


% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005/elecResp_n616_p45';
% elec = 43; %also try 42
% 
% elecRespPathBaseOtherSpikes = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005/elecResp_ndummy_p45';
% compareFreq = 10;
% comparisonResponse.movieInd = 9;
% nSequence = 20;


% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2010-10-28-3/data003/elecResp_n242_p15';
% elec = 17;
% 
% elecRespPathBaseOtherSpikes = '';
% compareFreq = 10;
% comparisonResponse.movieInd = 7;
% nSequence = 20;

% elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data004/elecResp_n753_p51';
% nSequence = 20;

elecRespPathBase = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data011/elecResp_n590_p40';
elec = 42;

elecRespPathBaseOtherSpikes = '';
compareFreq = 10;
comparisonResponse.movieInd = 7;
nSequence = 20;




% plotting parameters

rasterParams.rasterType = 'line'; %'dot' or 'line'

psthParams.markerSize = 5;

%time relative to stimulus where spikes where spikes were detected (in samples)
analyzedRegion = [8 35];



%% plot pulse-by-pulse response probabilities (over repetitions of pulse
% sequence) for a range of stimulus amplitudes

if 1
    movieInds = 10:15;
    
    %freq = 20;
    %elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
    %frequencyResponseHists(elecRespPath, freq, movieInds, nSequence, psthParams)
    
    %movieInds = 8:13;
    freq = 160;
    elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
    frequencyResponseHists(elecRespPath, freq, movieInds, nSequence, psthParams)
end

%% plot response probability vs. number of previous successes
if 1
    movieInds = 10:15;

    freq = 160;
    elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
    frequencyResponseHistsVsPrevSpikes(elecRespPath, freq, movieInds, nSequence, psthParams)
end


%% plot as rasters

if 0
    movieInd = 15;
    freq = 160;
    elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
    frequencyResponseRasters(elecRespPath, freq, movieInd, nSequence, analyzedRegion, rasterParams)
end

%% plot spike waveforms vs. repetition

%comparisonResponse.otherSpikePath = [elecRespPathBaseOtherSpikes '_f' num2str(compareFreq) '.mat'];
comparisonResponse.path = [elecRespPathBase '_f' num2str(compareFreq) '.mat'];
comparisonResponse.freq = compareFreq;

excludeOrSubtract = 'exclude';
plotVsLatency = false;
subtractPulseByPulseArtifact = true;
yPlotOffset = 150;
traceLength = 30;

nthPulseToCompare = 15:20; %mean uncontaminated spike waveforms corresponding to these pulses are plotted for comparison (0 means don't plot, -1 means plot comparison frequency)

if 1
    movieInd = 11;
    freq = 160;
    elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
    otherSpikePath = [elecRespPathBaseOtherSpikes '_f' num2str(freq) '.mat'];
    %otherSpikePath = '';
    
    frequencyResponseSpikeWaveforms(elecRespPath, freq, movieInd, nSequence, elec,...
        'otherSpikePath', otherSpikePath, 'comparisonResponse', comparisonResponse, 'excludeOrSubtract',...
        excludeOrSubtract, 'plotVsLatency', plotVsLatency, 'yPlotOffset', yPlotOffset, 'traceLengthPlot', traceLength, 'subtractPulseByPulseArtifact', subtractPulseByPulseArtifact, 'nthPulseToCompare', nthPulseToCompare)
end

%% plot heat map of latencies

freq = 160;
movieInd = 13;

elecRespPath = [elecRespPathBase '_f' num2str(freq) '.mat'];
latencies = extractFrequencyAnalysis(elecRespPath, movieInd, nSequence);

figure('position', [300 200 200 200])
imagesc(latencies', [5 25])
xlabel('pulse number')
ylabel('repetition')
title('latency heat map')


%% plot multiple response curves corresponding to different frequencies for
% same cell/stim electrode

if 0

frequencies = [5 10 20 40 80 160];
markerSize = 5;

nFreq = length(frequencies);
figure('position', [100 100 400 300])
hold on

plotColors = hsv(nFreq);
legendStrings = cell(1,nFreq);

for ii = 1:nFreq
    elecRespPathTemp = [elecRespPathBase '_f' num2str(frequencies(ii)) '.mat'];
    
    load(elecRespPathTemp);
    
    data = zeros(2, length(elecResp.stimInfo.movieNos));
    data(1,:) = elecResp.stimInfo.stimAmps;
    data(2,:) = elecResp.analysis.successRates;
    lockedAmps = elecResp.analysis.finalized;
    for jj = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{jj})
            data(:,jj) = [];
            lockedAmps(jj) = [];
        end
    end
    data(1,:) = abs(data(1,:));
    
    a = elecResp.analysis.erfParams(1);
    b = elecResp.analysis.erfParams(2);
    
    nAmps = size(data, 2);
    xProj = data(1,1):0.001:data(1,nAmps);
    projection = 0.5 + 0.5*erf(a*xProj+b);
    
    hold on
    
    for jj = 1:size(data, 2)
        if lockedAmps(jj)
            hTemp = plot(data(1,jj), data(2,jj),'.', 'markerSize', markerSize, 'markerEdgeColor', plotColors(ii,:));
        else
            hTemp = plot(data(1,jj), data(2,jj),'.', 'markerSize', markerSize, 'markerEdgeColor', [0.6 0.6 0.6]);
        end
        set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    end
    
    hTemp = plot(data(1,:), data(2,:), '--', 'color', plotColors(ii,:));
    set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    
    hTemp = plot([-b/a -b/a], [0 0.5], 'color', plotColors(ii,:));
    set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    
    
    plot(xProj, projection,'-', 'color', plotColors(ii,:));
    legendStrings{ii} = [num2str(frequencies(ii)) 'Hz'];
end

shortName = elecResp.names.rrs_short_name;
fIndex = strfind(shortName, '_f');
shortName = shortName(1:fIndex-1);
underIndex = strfind(shortName, '_');
for ii = 1:length(underIndex)
    shortName = [shortName(1:underIndex(ii)-1) ' ' shortName(underIndex(ii)+1:end)];
end

legend(legendStrings, 'location', 'SouthEast')

set(gca, 'ylim', [0 1])
xlabel('current amplitude (µA)')
ylabel('response rate')
title(['frequency-dependent responses ', 10, shortName])
hold off

end
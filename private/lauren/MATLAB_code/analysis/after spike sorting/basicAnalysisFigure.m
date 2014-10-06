clear all


%% off midget example

%pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data006/elecResp_n5_p1_50us';

%movieNumber = 52;


%% on midget example

%pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data006/elecResp_n31_p3_50us';

%movieNumber = 22; %40 is also good (19 has no evoked spikes)

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/elecResp_n858_p58';
movieNumber = 20;

%% on midget with good range of responses but eis not great toward end

%pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/elecResp_n782_p54';
%movieNumber = 13; %roughly 10% successes
% movieNumber = 24 %one failure

%% on parasol example that has good examples of 0%, 50% and 100% stimulation

%pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/elecResp_n257_p19';

%movieNumber = 4; % 0% response probability
%movieNumber = 17; % ~50% response probability
%movieNumber = 24; % 100% response probability



xRange = [0 2]; % in ms
%yRange = [-400 -50]; % for raw data
%yRangeSub = [-200 150]; %for subtracted data

yRange = [-400 -50]; % for raw data
yRangeSub = [-200 150]; %for subtracted data

%%
% gets data
temp = load(pathToData);
elecResp = temp.elecResp;

movieIndex = find(elecResp.stimInfo.movieNos == movieNumber);

dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
    elecResp.stimInfo.movieNos(movieIndex), 99999);

centerChannel = elecResp.cells.recElec;
nPulses = elecResp.stimInfo.nPulses(movieIndex);
latencies = [elecResp.analysis.latencies{movieIndex} elecResp.analysis.otherLatencies{movieIndex}];

% raw data
dataToPlot = squeeze(dataTraces(:, centerChannel, :));

%subtract estimated artifact
subtractionVector = elecResp.analysis.estArtifact{movieIndex}(elecResp.cells.goodElecs == centerChannel, :)';
elecResp.analysis.estArtifact{movieIndex}(elecResp.cells.goodElecs == elecResp.cells.recElec, :);

% after subtracting artifact
dataToPlotSubtracted = zeros(size(dataTraces, 1), size(dataTraces, 3));
for i = 1:nPulses
    dataToPlotSubtracted(i, :) = squeeze(dataTraces(i, centerChannel, :)) - subtractionVector;
end

neuronIDs = [elecResp.cells.main elecResp.cells.active{movieIndex}];
goodChannels = elecResp.cells.goodElecs;
centerChannelIndex = find(goodChannels == elecResp.cells.recElec);
nTemplates = length(neuronIDs);

figure
axes('position', [0.1 0.7 0.8 0.25])
%plotting successes
hold on
%plot failures in black
for i = 1:elecResp.stimInfo.nPulses(movieIndex)
    if ~any(latencies(i,:))
        plot(0.05:0.05:3.5, dataToPlot(i, :), 'k');
    end
end
%plot successes in red
for i = 1:elecResp.stimInfo.nPulses(movieIndex)
    if any(latencies(i,:))
        plot(0.05:0.05:3.5, dataToPlot(i, :), 'r')
    end
end
hold off
set(gca, 'xLim', xRange, 'xtick', [0 1 2 3], 'yLim', yRange)

axes('position', [0.1 0.4 0.8 0.25])
%plotting successes
hold on
%plot failures in black
for i = 1:elecResp.stimInfo.nPulses(movieIndex)
    if ~any(latencies(i,:))
        plot(0.05:0.05:3.5, dataToPlotSubtracted(i, :), 'k')
    end
end
%plot successes in red
for i = 1:elecResp.stimInfo.nPulses(movieIndex)
    if any(latencies(i,:))
        plot(0.05:0.05:3.5, dataToPlotSubtracted(i, :), 'r')
    end
end
hold off
set(gca, 'xLim', xRange, 'xtick', [0 1 2 3], 'yLim', yRangeSub)


% raster/psth plot
rasterOptions.markerSize = 5; %for raster plot
rasterOptions.markerColor = [0 0 0];
psthOptions.binSize = 0.05; %in ms
psthOptions.lineColor = [0 0 0];

figure('position', [100 100 200 200])
axes('units', 'pixels', 'position', [20 20 161.8 100])
plot_raster(gca, elecResp, movieNumber, rasterOptions)
xlabel('')
set(gca, 'xLim', xRange, 'xtick', [])

figure('position', [100 100 200 200])
axes('units', 'pixels', 'position', [20 20 161.8 100])
plot_psth(gca, elecResp, movieNumber, psthOptions)
set(gca, 'xLim', xRange)




%% from refreshEiPlots

break

templateColors = hsv(nTemplates);

showTemplates = 1;


eiData = cell(nTemplates, 1); %stores EIs of active neurons: eiData{neuronIndex}(channels, samples)

eiData{1} = elecResp.cells.mainEI;
for i = 2:nTemplates
    eiData{i} = elecResp.cells.allEIs{elecResp.cells.all == elecResp.cells.active{movieIndex}(i-1)};
end

eiYMin = 5000; eiYMax = -5000;
for i = 1:nTemplates
    eiYMin = min([min(min(eiData{i})), eiYMin]);
    eiYMax = max([max(max(eiData{i})), eiYMax]);
end

failuresBin = (latencies == 0);
failuresBin = min(failuresBin, [], 2);

if sum(failuresBin) %some traces have no spikes
    channelsToUse = getCluster(centerChannel);
    %         for i = length(channelsToUse):-1:1
    %             if channelsToUse(i) == 4
    %                 channelsToUse(i) = [];
    %             end
    %         end
else %all traces have spikes
    channelsToUse = elecResp.cells.goodElecs;
end

nChannels = length(channelsToUse);


ei = cell(nTemplates, 1);
targetEi = cell(nTemplates, 1);
targetEiMinPos = zeros(nTemplates, 1);
nSpikes = zeros(nTemplates, 1);
for i = 1:nTemplates
    ei{i} = findStimEi(dataTraces, elecResp, movieNumber, i-1, eiData, channelsToUse, false);
    targetEi{i} = eiData{i}(channelsToUse, :);
    targetEiMinPos(i) = find(squeeze(targetEi{i}(channelsToUse==centerChannel,:))...
        ==min(squeeze(targetEi{i}(channelsToUse==centerChannel,:)))); %position of minimum on primary electrode
    nSpikes(i) = sum(latencies(:,i)~=0);

    eiYMin = min([min(min(min(ei{i}))), eiYMin]);
    eiYMax = max([max(max(max(ei{i}))), eiYMax]);
end


figure('position', [100 100 800 200])

for j = 1:7
    axes('units', 'pixels', 'position', [50 + 100*(j-1), 30, 90, 150])
    hold on

    elecIndex = j;

    for i = 1:nTemplates
        if nSpikes(i)
            for k = 1:nPulses
                if latencies(k,i)
                    plot((1:26)/20, squeeze(ei{i}(k,elecIndex,:)), 'Color', 'r');
                end
            end
        end
        if showTemplates
            current = plot((1:26)/20, squeeze(targetEi{i}(elecIndex, targetEiMinPos(i)-10:targetEiMinPos(i)+15)), 'LineWidth', 2, 'Color', 'k');
        end
    end
    hold off
    set(gca, 'ylim', [eiYMin eiYMax], 'xlim', [0 1.5], 'ytick', [])
end













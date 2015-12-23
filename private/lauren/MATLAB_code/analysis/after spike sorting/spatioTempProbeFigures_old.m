%script to make summary plots from spatiotemporal probe data

%preElecs = [16 13 11 12 17 18];
preElecs = [13 11 12 17 18];
neuronID = 999;
stimElec = 14;
offsets = [0 10 20 40 80 160];
pathToElecResp = '/snle/lab/Experiments/Array/Analysis/2010-03-05-4/data006/';
latencyAmps = [-0.55606 -0.60661 -0.68244]; %amplitudes at which to measure latency

temp = load([pathToElecResp 'eiFile998.mat']);
eiData{1} = zeros(64, 1);
for i = 1:64
    if ~(i==9||i==25||i==57)
        eiData{1}(i) = max(abs(temp.ei(i,:))); %#ok<AGROW>
    end
end
eiData{1}(31) = 0;

% plotting options
threshYLims = [0.4 0.7];

binSize = 0.1; % in ms
binEnd = 1.5; %in ms

lineColors(1,:) = [90 156 0]/255; %pale grass
lineColors(2,:) = [255 124 59]/255; %salmon
lineColors(3,:) = [101 52 255]/255; %purple
lineColors(4,:) = [52 198 247]/255; %aqua
lineColors(5,:) = [238 55 128]/255; %calm magenta


nPreElecs = length(preElecs);
nOffsets = length(offsets);
nLatencies = length(latencyAmps);

binEdges = 0:binSize:binEnd;

%% special case: some successes to prepulse

latenciesAfterFail = cell(1, nOffsets, nLatencies);
latenciesAfterSucc = cell(1, nOffsets, nLatencies);
thresholdsAfterFail = zeros(1, nOffsets);
thresholdsAfterSucc = zeros(1, nOffsets);
threshStdAfterFail = zeros(1, nOffsets);
threshStdAfterSucc = zeros(1, nOffsets);


for j = 3:nOffsets
    stimPulsePath = [pathToElecResp '/elecResp_n999_p14_pre16_d' num2str(offsets(j)) '.mat'];
    temp = load([pathToElecResp '/elecResp_n999_p14_pre16_d' num2str(offsets(j)) '_preResponse.mat']);
    elecRespPrePulse = temp.elecResp;
    movieIndeces = zeros(nLatencies, 1);
    
    
    nMovies = length(elecRespPrePulse.stimInfo.movieNos);
    subsetBinFail = cell(nMovies, 1);
    subsetBinSucc = cell(nMovies, 1);
    for i = 1:nMovies
        if ~isempty(elecRespPrePulse.analysis.type{i})
            subsetBinFail{i} = elecRespPrePulse.analysis.latencies{i}==0; %failures to prepulse
            subsetBinSucc{i} = elecRespPrePulse.analysis.latencies{i}~=0; %successes to prepulse
        end
    end
    [latenciesAfterFailTemp thresholdsAfterFail(j) threshStdAfterFail(j)] = fitErfToSubset(stimPulsePath, subsetBinFail);
    [latenciesAfterSuccTemp thresholdsAfterSucc(j) threshStdAfterSucc(j)] = fitErfToSubset(stimPulsePath, subsetBinSucc);
    
    
    temp = load([pathToElecResp 'elecResp_n999_p14_pre16_d' num2str(offsets(j))]);
    elecResp = temp.elecResp;
    
    for k = 1:nLatencies
        movieIndex = find((abs(elecResp.stimInfo.stimAmps) > abs(0.99*latencyAmps(k))) .* (abs(elecResp.stimInfo.stimAmps) < abs(1.01*latencyAmps(k))));
        latenciesAfterFail{1, j, k} = latenciesAfterFailTemp{movieIndex};
        latenciesAfterSucc{1, j, k} = latenciesAfterSuccTemp{movieIndex};
    end
end

%% extract data: normal cases

maxNPulses = 0;

thresholds = zeros(nPreElecs, nOffsets);
thresholdStds = zeros(nPreElecs, nOffsets);

latencies = cell(nPreElecs, nOffsets, nLatencies);

for i = 1:nPreElecs
    for j = 1:nOffsets
        temp = load([pathToElecResp 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec)...
            '_pre' num2str(preElecs(i)) '_d' num2str(offsets(j))]);
        elecResp = temp.elecResp;
        
        thresholds(i,j) = elecResp.analysis.threshold;
        thresholdStds(i,j) = elecResp.analysis.threshStd;
        
        for k = 1:nLatencies
            movieIndex = find((abs(elecResp.stimInfo.stimAmps) > abs(0.99*latencyAmps(k))) .* (abs(elecResp.stimInfo.stimAmps) < abs(1.01*latencyAmps(k))));
            latencies{i,j,k} = elecResp.analysis.latencies{movieIndex};
        end
    end
end

% stim pulse without prepulse
temp = load([pathToElecResp 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec)]);
elecResp = temp.elecResp;
thresholdIsolated = elecResp.analysis.threshold;
thresholdStdIsolated = elecResp.analysis.threshStd;

latenciesIsolated = cell(nLatencies, 1);
for k = 1:nLatencies
    movieIndex = find((abs(elecResp.stimInfo.stimAmps) > abs(0.99*latencyAmps(k))) .* (abs(elecResp.stimInfo.stimAmps) < abs(1.01*latencyAmps(k))));
    latenciesIsolated{k} = elecResp.analysis.latencies{movieIndex};
end

    
%% plotting: special case

figure('position', [100 100 1000 600], 'color', [1 1 1])

%ei plot
eiAxes = axes('units', 'pixels', 'position', [10 410 180 180]);
plotEi61('', neuronID, 'axesH', eiAxes, 'eiData', eiData, 'markElecs', [stimElec 16])

%thresholds plot
threshAxes = axes('units', 'pixels', 'position', [300 450 300 140]);
hold on
plot(offsets(3:end)/20, thresholdsAfterFail(1, 3:end), 'Color', lineColors(4,:))
plot(offsets(4:end)/20, thresholdsAfterSucc(1, 4:end), 'Color', lineColors(5,:))
for j = 3:nOffsets
    plot([offsets(j)/20 offsets(j)/20],...
        [thresholdsAfterFail(1,j) - threshStdAfterFail(1,j) thresholdsAfterFail(1,j) + threshStdAfterFail(1,j)],...
        'Color', lineColors(4,:))
    if j>3
        plot([offsets(j)/20 offsets(j)/20],...
            [thresholdsAfterSucc(1,j) - threshStdAfterSucc(1,j) thresholdsAfterSucc(1,j) + threshStdAfterSucc(1,j)],...
            'Color', lineColors(5,:))
    end
end
plot(offsets(end)/20 + 5, thresholdIsolated, 'k.', 'markerSize', 3)
plot([offsets(end)/20 + 5 offsets(end)/20 + 5], [thresholdIsolated - thresholdStdIsolated thresholdIsolated + thresholdStdIsolated], 'k')

set(threshAxes, 'xlim', [offsets(1)/20-0.5 offsets(end)/20+5.5], 'ylim', threshYLims, 'xtick', offsets(1)/20:2:offsets(end)/20)
hold off
ylabel('threshold (µA)')
xlabel('delay between prepulse and stim pulse (ms)')

%legends
legendAxes = axes('units', 'pixels', 'position', [700 450 200 150]);
hold on
plot(0.1, 0.6, 's', 'MarkerSize', 5,...
    'MarkerFaceColor', lineColors(4,:), 'MarkerEdgeColor', lineColors(4,:))
plot(0.1, 0.4, 's', 'MarkerSize', 5,...
    'MarkerFaceColor', lineColors(5,:), 'MarkerEdgeColor', lineColors(5,:))
text(0.2, 0.6, 'after no response to prepulse')
text(0.2, 0.4, 'after response to prepulse')

set(legendAxes, 'ylim', [0 1], 'xlim', [0 1])
hold off
axis off

legendAxes2 = axes('units', 'pixels', 'position', [50 250 200 100]);
hold on
for k = 1:nLatencies
    plot(0.1, (nLatencies-k)*(0.8/nLatencies) + 0.1, 's', 'MarkerSize', 5,...
        'MarkerFaceColor', lineColors(k,:), 'MarkerEdgeColor', lineColors(k,:))
    text(0.2, (nLatencies-k)*(0.8/nLatencies) + 0.1, [num2str(latencyAmps(k)) ' µA'])
end
set(legendAxes, 'ylim', [0 1], 'xlim', [0 1])
hold off
axis off

%latency histograms
lAxes = cell(nOffsets+1, 2);
for j = 3:nOffsets
    lAxes{j,1} = axes('units', 'pixels', 'position', [50+(j-1)*900/(nOffsets+1) 250 800/(nOffsets+1) 120]);
    title([num2str(offsets(j)/20) 'ms delay'])
    lAxes{j,2} = axes('units', 'pixels', 'position', [50+(j-1)*900/(nOffsets+1) 50 800/(nOffsets+1) 120]);
    hold on
    for k = 1:nLatencies
        if ~isempty(latenciesAfterFail{1,j,k})
            psthPlotterBase(lAxes{j,1}, binEdges, (latenciesAfterFail{1,j,k} - offsets(j))/20,...
                'lineColor', lineColors(k,:), 'fillHist', true, 'normalize', true)
        end

        if ~isempty(latenciesAfterSucc{1,j,k})
            psthPlotterBase(lAxes{j,2}, binEdges, (latenciesAfterSucc{1,j,k} - offsets(j))/20,...
                'lineColor', lineColors(k,:), 'fillHist', true, 'normalize', true)
        end
    end
    hold off
    set(lAxes{j,1}, 'ylim', [0 1])
    set(lAxes{j,2}, 'ylim', [0 1])
end
lAxes{nOffsets+1, 1} = axes('units', 'pixels', 'position', [50+(nOffsets)*900/(nOffsets+1) 150 800/(nOffsets+1) 120]);

hold on
for k = 1:nLatencies
    if ~isempty(latenciesIsolated{k})
        psthPlotterBase(lAxes{nOffsets+1, 1}, binEdges, (latenciesIsolated{k})/20,...
            'lineColor', lineColors(k,:), 'fillHist', true, 'normalize', true)
    end
end
hold off
set(lAxes{nOffsets+1, 1}, 'ylim', [0 1])
title('no prepulse')

xlabel(lAxes{3,2}, 'time (ms)')
ylabel(lAxes{3,2}, 'spike probability')
axes(lAxes{3,1})
text(0.2, 0.9, ['after no response' 10 'to prepulse'])
axes(lAxes{3,2})
text(0.2, 0.9, ['after response' 10 'to prepulse'])


%% plot

for i = 1:nPreElecs
    figure('position', [100 100 1000 400], 'color', [1 1 1], 'name', ['prepulse on electrode' num2str(preElecs(i))])
    
    %ei plot
    eiAxes = axes('units', 'pixels', 'position', [10 210 180 180]);
    plotEi61('', neuronID, 'axesH', eiAxes, 'eiData', eiData, 'markElecs', [stimElec preElecs(i)])
    
    %thresholds plot
    threshAxes = axes('units', 'pixels', 'position', [300 250 300 140]);
    hold on
    plot(offsets/20, thresholds(i, :), 'k')
    plot(offsets(end)/20 + 5, thresholdIsolated, 'k.', 'markerSize', 3)
    set(threshAxes, 'xlim', [offsets(1)/20-0.5 offsets(end)/20+5.5], 'ylim', threshYLims, 'xtick', offsets(1)/20:2:offsets(end)/20)
    for j = 1:nOffsets
        plot([offsets(j)/20 offsets(j)/20],...
            [thresholds(i,j) - thresholdStds(i,j) thresholds(i,j) + thresholdStds(i,j)], 'k')
    end
    plot([offsets(end)/20 + 5 offsets(end)/20 + 5], [thresholdIsolated - thresholdStdIsolated thresholdIsolated + thresholdStdIsolated], 'k')
    hold off
    ylabel('threshold (µA)')
    xlabel('delay between prepulse and stim pulse (ms)')
    
    %legend
    legendAxes = axes('units', 'pixels', 'position', [700 250 200 150]);
    hold on
    for k = 1:nLatencies
        plot(0.1, (nLatencies-k)*(0.8/nLatencies) + 0.1, 's', 'MarkerSize', 5,...
            'MarkerFaceColor', lineColors(k,:), 'MarkerEdgeColor', lineColors(k,:))
        text(0.2, (nLatencies-k)*(0.8/nLatencies) + 0.1, [num2str(latencyAmps(k)) ' µA'])
    end
    set(legendAxes, 'ylim', [0 1], 'xlim', [0 1])
    hold off
    axis off
    
    %latency histograms
    lAxes = cell(nOffsets+1, 1);
    for j = 1:nOffsets
        lAxes{j} = axes('units', 'pixels', 'position', [50+(j-1)*900/(nOffsets+1) 50 800/(nOffsets+1) 120]);
        hold on
        for k = 1:nLatencies
            if ~isempty(latencies{i,j,k})
                psthPlotterBase(lAxes{j}, binEdges, (latencies{i,j,k} - offsets(j))/20,...
                    'lineColor', lineColors(k,:), 'fillHist', true, 'normalize', true)
            end
        end
        hold off            
        set(lAxes{j}, 'ylim', [0 1])
        title([num2str(offsets(j)/20) 'ms delay'])
    end    
    xlabel(lAxes{1}, 'time (ms)')
    ylabel(lAxes{1}, 'spike probability')
    
    
    lAxes{nOffsets+1} = axes('units', 'pixels', 'position', [50+(nOffsets)*900/(nOffsets+1) 50 800/(nOffsets+1) 120]);

    hold on
    for k = 1:nLatencies
        if ~isempty(latenciesIsolated{k})
            psthPlotterBase(lAxes{nOffsets+1, 1}, binEdges, (latenciesIsolated{k})/20,...
                'lineColor', lineColors(k,:), 'fillHist', true, 'normalize', true)
        end
    end
    hold off
    set(lAxes{nOffsets+1}, 'ylim', [0 1])
    title('no prepulse')
    
    
end



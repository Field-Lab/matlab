function clusterIndeces = manualClustering()


%% script to run hand-clustering

clear all
close all

DataPath = '/Analysis.noindex/Lauren/2008-08-26-0/data006_proba';
nPatterns = 1;
nMovies = 1;

dataTraces = cell(nPatterns, nMovies);
channels = cell(nPatterns, nMovies);

PCAWindowStart = 5;
PCAWindowEnd = 30;



%% Construct the gui components

h.gui = figure('Visible', 'off', 'Position',[100 100 800 600], 'Name', 'Hand-Clustering Main.', 'Resize', 'off');

h.p1h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [50 300 200 200]); %from left, from bottom, width, height
h.p2h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [300 300 200 200]);
h.p3h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [550 300 200 200]);
h.p4h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [50 50 200 200]);
h.p5h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [300 50 200 200]);
h.p6h = axes('Parent', h.gui, 'Units', 'pixels', 'Position', [550 50 200 200]);


uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC2', 'Position', [50 500 200 30], 'Callback', @PC12);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC3', 'Position', [300 500 200 30], 'Callback', @PC13);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC1 vs PC4', 'Position', [550 500 200 30], 'Callback', @PC14);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC2 vs PC3', 'Position', [50 250 200 30], 'Callback', @PC23);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC2 vs PC4', 'Position', [300 250 200 30], 'Callback', @PC24);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'PC3 vs PC4', 'Position', [550 250 200 30], 'Callback', @PC34);
uicontrol(h.gui, 'Style', 'pushbutton', 'String', 'cluster all remaining', 'Position', [300 550 200 30], 'Callback', @clusterAll);

%% initialization tasks

axes(h.p1h)
plot(PCAScore(:,1), PCAScore(:,2),'k.')
% xlabel('PC1')
% ylabel('PC2')

%% Callbacks for guih

function PC12(hObject, eventdata)
    PCx = 1;
    PCy = 2;
    close;
end



%% loading data
for i = 1:nPatterns
    if i == 9 || i == 25 || i == 57 || i == 4
    else
        disp(sprintf('loading pattern %.0f',i))
        for j = 1:nMovies
            [dataTracesFull, channelsFull] = NS_ReadPreprocessedDataNoArtifacts(DataPath,i,j);
            electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
            surroundingChannels = electrodeMap.getAdjacentsTo(i, 1)';
            surroundingChannels = surroundingChannels(surroundingChannels ~= 4); %removes channel 4
            dataTraces{i,j} = dataTracesFull(:,surroundingChannels,:);
            channels{i,j} = channelsFull(surroundingChannels);
        end
    end
end

clusterIndex = cell(nPatterns, nMovies);

i = 1;
j = 1;
clusterIndex{i,j} = manualPCACluster(dataTraces{i,j}(:,:,PCAWindowStart:PCAWindowEnd));


% plots representation of clusterIndex
clusterColor = hsv(max(clusterIndex{i,j}));
nTraces = size(dataTraces{i,j}, 1);
nElectrodes = size(dataTraces{i,j}, 2);
nPCASamples = PCAWindowEnd - PCAWindowStart + 1;


PCAData = dataTraces{i,j}(:,:,PCAWindowStart:PCAWindowEnd);
prinCompArray = zeros(nTraces, nElectrodes*nPCASamples);
for k = 1:nElectrodes % concatenates traces on different electrodes electrode for each pulse
    prinCompArray(:, (k-1)*nPCASamples + 1 : k*nPCASamples) = PCAData(:,k,:);
end

[PCACoef, PCAScore] = princomp(prinCompArray);


traceIndeces = lassoTraces(prinCompArray);

selectedTraces = zeros(1, nTraces);
selectedTraces(traceIndeces) = 1;

figure %clusters based on lassoTraces

% plots traces as different colors representing the cluster they belong to.
for k = 1:length(clusterIndex{i,j})
    subplot(2,1,1)
    hold on
    if selectedTraces(k)
        plot(PCAScore(k,1), PCAScore(k,2), 'm.', 'MarkerSize', 20)
    else
        plot(PCAScore(k,1), PCAScore(k,2), 'k.', 'MarkerSize', 20)
    end
    hold off
    subplot(2,1,2)
    hold on
    if selectedTraces(k)
        plot(prinCompArray(k,:), 'm-');
    else
        plot(prinCompArray(k,:), 'k-');
    end
    hold off
end


% figure
%
% % plots traces as different colors representing the cluster they belong to.
% for k = 1:length(clusterIndex{i,j})
%     subplot(2,1,1)
%     hold on
%     plot(PCAScore(k,1), PCAScore(k,2), '.', 'MarkerFaceColor', clusterColor(clusterIndex{i,j}(k),:), 'MarkerEdgeColor', clusterColor(clusterIndex{i,j}(k),:), 'MarkerSize', 20)
%     hold off
%     subplot(2,1,2)
%     hold on
%     current = plot(prinCompArray(k,:));
%     set(findobj(current,'Type','line'), 'Color', clusterColor(clusterIndex{i,j}(k),:))
%     hold off
% end

end %end of main function

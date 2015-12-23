clear all

dataPath = '/snle/lab/Experiments/Array/Analysis/2011-10-24/data002/';
arrayName = 'KM43';
platType = 'CC+CV';

patternNos = [1:8 10:24 26:56 58:64];


currentMovie = 1;

stimSignal = zeros(64, 100);
for ii = 1:length(patternNos)
    dataTraces=NS_ReadPreprocessedData([dataPath filesep 'p' num2str(patternNos(ii))], '', 0, patternNos(ii), currentMovie, 99999);
    signalTemp = squeeze(mean(dataTraces,1));
    stimSignal(patternNos(ii),:) = signalTemp(patternNos(ii),:);
end


yMax = max(max(stimSignal(stimSignal ~= 0)));
yMin = min(min(stimSignal(stimSignal ~= 0)));

maxRange = 0;
ranges = zeros(1,64);
for ii = 1:64
    ranges(ii) = max(stimSignal(ii,:))-min(stimSignal(ii,:));
    maxRange = max(maxRange, ranges(ii));
end

rangesNorm = ranges/maxRange;

[xCoords yCoords] = getElectrodeCoords61();


figure('position', [100 100 600 525])
hold on
axis equal

for i = 1:64
    if ~(i==9||i==25||i==57)
        if rangesNorm(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(rangesNorm(i)*40), 'MarkerFaceColor', 'k')
        end
    end
end
axis off

axes('position', [0.05 0.9 0.3 0.1])
text(0, 0.5, [arrayName, 10, 'platinized with ' platType, 10 'max signal range = ' num2str(maxRange) 10 'stim elec artifacts'])
axis off

%%

neighborSignalSum = cell(64,1);
for ii = 1:64
    neighborSignalSum{ii} = zeros(1,100);
end

neighborSignalN = zeros(64,1);

for ii = 1:length(patternNos)
    
    patternNo = patternNos(ii);
    
    
    dataTraces=NS_ReadPreprocessedData([dataPath filesep 'p' num2str(patternNo)], '', 0, patternNo, currentMovie, 99999);
    
    signalTemp = squeeze(mean(dataTraces,1));
    
    %store mean signal on electrodes neighboring stim electrode
    neighbors = getCluster(patternNo);
    neighbors = neighbors (2:end);
    
    %keyboard
    for jj = 1:length(neighbors)
        if ~isempty(neighborSignalSum{jj})
            neighborSignalSum{neighbors(jj)} = neighborSignalSum{neighbors(jj)} + signalTemp(neighbors(jj),:);
        else
            neighborSignalSum{neighbors(jj)} = signalTemp(neighbors(jj),:);
        end
        neighborSignalN(neighbors(jj)) = neighborSignalN(neighbors(jj)) + 1;
    end    
end

    
signal = zeros(64, 100);
for ii = 1:64
    signal(ii,:) = neighborSignalSum{ii}/neighborSignalN(ii);
end
    
%% generate plots

[xCoords yCoords] = getElectrodeCoords61();

arrayWidth = max(xCoords)*2*1.2;
arrayHeight = max(yCoords)*2*1.2;

figure('position', [100 100 600 525])

yMax = max(max(signal(signal ~= 0)));
yMin = min(min(signal(signal ~= 0)));

maxRange = 0;
ranges = zeros(1,64);
for ii = 1:64
    ranges(ii) = max(signal(ii,:))-min(signal(ii,:));
    maxRange = max(maxRange, ranges(ii));
end

rangesNorm = ranges/maxRange;

for i = 1:64
    if ~(i==9||i==25||i==57)
        axes('position', [xCoords(i)/arrayWidth + 0.45, yCoords(i)/arrayHeight + 0.45, 0.1, 0.1])
        plot(signal(i,1:20), 'k')
        set(gca, 'ylim', [yMin yMax])
        %set(gca, 'ylim', [-10 8], 'xtick', [], 'ytick', [])
        axis off
    end
end

axes('position', [0.05 0.9 0.3 0.1])
text(0, 0.5, [arrayName, 10, 'platinized with ' platType, 10 'max signal range = ' num2str(maxRange)])
axis off

% array orientation legend
axes('position', [0.8 0.8 0.15 0.15])
hold on
plot(xCoords([8 19 30 40 51 62 8]), yCoords([8 19 30 40 51 62 8]), 'k-')
text(xCoords(7),  yCoords(7),  '8',  'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(61), yCoords(61), '62', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(50), yCoords(50), '51', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(39), yCoords(39), '40', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(29), yCoords(29), '30', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(18), yCoords(18), '19', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
hold off
axis equal
axis off




%% figure using dots only

figure('position', [100 100 600 525])
hold on
axis equal

for i = 1:64
    if ~(i==9||i==25||i==57)
        if rangesNorm(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(rangesNorm(i)*40), 'MarkerFaceColor', 'k')
        end
    end
end
axis off

axes('position', [0.05 0.9 0.3 0.1])
text(0, 0.5, [arrayName, 10, 'platinized with ' platType, 10 'max signal range = ' num2str(maxRange)])
axis off


% array orientation legend
axes('position', [0.8 0.8 0.15 0.15])
hold on
plot(xCoords([8 19 30 40 51 62 8]), yCoords([8 19 30 40 51 62 8]), 'k-')
text(xCoords(7),  yCoords(7),  '8',  'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(61), yCoords(61), '62', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(50), yCoords(50), '51', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(39), yCoords(39), '40', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(29), yCoords(29), '30', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
text(xCoords(18), yCoords(18), '19', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle')
hold off
axis equal
axis off
    
    
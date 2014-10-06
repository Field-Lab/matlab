%function figHandle = summaryPlotGenerator(dataTraces, centerChannel, channelRadius, data,...
%    artifact, template, results, shift)

function figHandle = summaryPlotGenerator(elecResp, dataTraces, channelRadius, movieNumber,...
    artifact, templates, results)

centerChannel = elecResp.cells.recElec;
goodChannels = elecResp.cells.goodElecs;
dataPath = elecResp.names.data_path;
patternNumber = elecResp.stimInfo.patternNo;
modelType = elecResp.analysis.type{elecResp.stimInfo.movieNos == movieNumber};
pathToEi = elecResp.names.rrs_ei_path;
neuronIDs = elecResp.cells.active;

actualArtifact = artifact.actual;
actualArtifactInit = artifact.actualInit;
artifactMean = artifact.mean;
if isfield(artifact, 'estimate')
    estArtifact = artifact.estimate;
end

latenciesFull = results.latenciesFull;
successesInit = results.successesInit;
successesFull = results.successesFull;


nTraces  = size(dataTraces, 1);
%nSamples = size(dataTraces, 3);
nTemplates = length(templates);
nElecs = length(goodChannels);

artFromSubset = false; %whether (false) or not (true) to include failure traces that have analysis flags when calculating the estimated artifact


%% figuring out good axis limits

ymin1 = min(min(squeeze(dataTraces(:,centerChannel,:))));
ymax1 = max(max(squeeze(dataTraces(:,centerChannel,:))));
ymin1 = 50*floor(ymin1/50);
ymax1 = 50*ceil(ymax1/50);

if ~any(any(isnan(actualArtifactInit))) && ~any(any(isnan(actualArtifact)))
    ymin2 = min([min(min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifactInit))...
        min(min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifact))]);
    ymax2 = max([max(max(reshape(max(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifactInit))...
        max(max(reshape(max(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifact))]);
elseif ~any(any(isnan(actualArtifact)))
    ymin2 = min(min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifact));
    ymax2 = max(max(reshape(max(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifact));
elseif ~any(any(isnan(actualArtifactInit)))
    ymin2 = min(min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifactInit));
    ymax2 = max(max(reshape(max(dataTraces(:, goodChannels, :)), length(goodChannels), []) - actualArtifactInit));
else
    ymin2 = min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []));
    ymax2 = min(reshape(min(dataTraces(:, goodChannels, :)), length(goodChannels), []));
end
ymin2 = 50*floor(ymin2/50);
ymax2 = 50*ceil(ymax2/50);


xPlotRange = 40;

figHandle = figure('Position',[100, 100, 800, 800]);

%% initial portion of algorithm

%%%%%%%%%%%%%% plots traces with spikes
subplot(6,2,1)

stimAmp = max(abs(getStimAmps(dataPath, patternNumber, movieNumber)));
text(-0.2, 1.5, [dataPath, 10,...
    'pattern ', num2str(patternNumber), ', movie ',...
	num2str(movieNumber), ', current amp: ', num2str(stimAmp),...
    ', ', modelType], 'Units','normalized')

hold on
nSuccesses = 0;
if sum(reshape(successesInit,[],1))
    successesToPlot = [];
    failuresToPlot = [];
    for i = 1:nTraces
        if any(successesInit(i,:))
            successToPlot = [];
            for j = 1:nElecs
                successToPlot = [successToPlot; squeeze(dataTraces(i,goodChannels(j),1:xPlotRange))]; %#ok<AGROW>
            end
            successesToPlot = [successesToPlot successToPlot]; %#ok<AGROW>
            nSuccesses = nSuccesses + 1;
        else
            failureToPlot = [];
            for j = 1:nElecs
                failureToPlot = [failureToPlot; squeeze(dataTraces(i,goodChannels(j),1:xPlotRange))]; %#ok<AGROW>
            end
            failuresToPlot = [failuresToPlot failureToPlot]; %#ok<AGROW>
        end
    end
    if ~isempty(failuresToPlot)
        current = plot(failuresToPlot);
        set(findobj(current,'Type','line'),'Color',[0.8 0.8 0.8])
    end
    plot(successesToPlot)
end

if exist('estArtifact', 'var') && max(max(abs(estArtifact)))
    estArtToPlot = [];
    for j = 1:nElecs
        estArtToPlot = [estArtToPlot estArtifact(j,1:xPlotRange)]; %#ok<AGROW>
    end
    plot(estArtToPlot,'k','LineWidth',2)
end
hold off
title('successes: raw data (partial algorithm)')
set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin1 ymax1])

%%%%%%%%%%%%%% plots traces with no spikes
subplot(6,2,3)
if nSuccesses ~= nTraces
    plot(failuresToPlot)
end
title('failures: raw data')
set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin1 ymax1])

%%%%%%%%%%%%%% plots traces with spikes after artifact subtraction
subplot(6,2,5)
if sum(reshape(successesInit,[],1)) && sum(reshape(~successesInit,[],1)) %not 0% or 100% success rate
    hold on
    for i = 1:nTraces
        if max(successesInit(i,:))
            spikesOnly = reshape(dataTraces(i, goodChannels, :), length(goodChannels), []) - actualArtifactInit;
            spikesOnlyCat = [];
            for j = 1:nElecs
                spikesOnlyCat = [spikesOnlyCat spikesOnly(j,1:xPlotRange)]; %#ok<AGROW>
            end
            plot(spikesOnlyCat, 'm')
        end
    end

    hold off
    title(['successes after subtraction of mean failure trace (' num2str(nSuccesses) ')'])
    set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin2 ymax2])
    xlabel('samples')
end



%% full algorithm

%%%%%%%%%%%%%% plots traces with spikes
subplot(6,2,2)
hold on

failuresBin = zeros(1, nTraces);
nSuccesses = 0;


if sum(reshape(successesFull,[],1))
    successesToPlotFull = [];
    failuresToPlotFull = [];
    for i = 1:nTraces
        if any(successesFull(i,:))
            successToPlotFull = [];
            for j = 1:nElecs
                successToPlotFull = [successToPlotFull; squeeze(dataTraces(i,goodChannels(j),1:xPlotRange))]; %#ok<AGROW>
            end
            successesToPlotFull = [successesToPlotFull successToPlotFull]; %#ok<AGROW>            
            nSuccesses = nSuccesses + 1;
             
        else
            failureToPlotFull = [];
            for j = 1:nElecs
                failureToPlotFull = [failureToPlotFull; squeeze(dataTraces(i,goodChannels(j),1:xPlotRange))]; %#ok<AGROW>
            end
            failuresToPlotFull = [failuresToPlotFull failureToPlotFull]; %#ok<AGROW>
            failuresBin(i) = 1;
            
        end
    end
    if ~isempty(failuresToPlotFull)
        current = plot(failuresToPlotFull);
        set(findobj(current,'Type','line'),'Color',[0.8 0.8 0.8])
    end
    plot(successesToPlotFull)
end

artMeanToPlot = [];
for j = 1:nElecs
    artMeanToPlot = [artMeanToPlot artifactMean(j,1:xPlotRange)]; %#ok<AGROW>
end
plot(artMeanToPlot,'k','LineWidth',2)

hold off

title('successes: raw data (full algorithm)')
set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin1 ymax1])


%%%%%%%%%%%%%% plots traces with no spikes
subplot(6,2,4)
if nSuccesses ~= nTraces
    plot(failuresToPlotFull)
end
title('failures: raw data')
set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin1 ymax1])



%%%%%%%%%%%%%% plots traces with spikes after artifact subtraction
subplot(6,2,6)
if sum(reshape(successesFull,[],1)) && sum(reshape(~successesFull,[],1)) %not 0% or 100% success rate
    hold on
    for i = 1:nTraces
        if max(successesFull(i,:))
            spikesOnly = reshape(dataTraces(i, goodChannels, :), length(goodChannels), []) - actualArtifactInit;
            spikesOnlyCat = [];
            for j = 1:nElecs
                spikesOnlyCat = [spikesOnlyCat spikesOnly(j,1:xPlotRange)]; %#ok<AGROW>
            end
            plot(spikesOnlyCat, 'r')
        end
    end
    hold off
    title(['successes after subtraction of mean failure trace (' num2str(nSuccesses) ')'])
    set(gca, 'xlim', [1 xPlotRange*nElecs], 'ylim', [ymin2 ymax2])
    xlabel('samples')
end


%% plots subtracted ei on surrounding


eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
eiData = cell(nTemplates, 1); %stores EIs of active neurons: eiData{neuronIndex}(channels, samples)
for i = 1:nTemplates
    eiData{i} = eiFile.getImage(neuronIDs(i));
    eiData{i} = squeeze(eiData{i}(1, 2:end, :));
end


failuresBin = (latenciesFull == 0);
failuresBin = min(failuresBin, [], 2);

if sum(failuresBin) %some traces have no spikes
    electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
    channelsToUse=electrodeMap.getAdjacentsTo(centerChannel, channelRadius)';

    % get rid of channel 4 (broken)
    for i = 1:length(channelsToUse)
        if channelsToUse(i) == 4;
            channelsToUse(i) = [];
            break
        end
    end
else%all traces have spikes
    channelsToUse = goodChannels;
end
nChannels = length(channelsToUse);


for i = 1:nTemplates
    
    ei = findStimEi(dataTraces, elecResp, movieNumber, i, eiData, channelsToUse, artFromSubset);
    
    %targetEi = reshape(eiData(1, channelsToUse+1, :), length(channelsToUse), []);
    targetEi = eiData{i}(channelsToUse, :);
    targetEiMinPos = find(squeeze(targetEi(channelsToUse==centerChannel,:))...
        ==min(squeeze(targetEi(channelsToUse==centerChannel,:)))); %position of minimum on primary electrode



    clear eiFile
    nSpikes = sum(latenciesFull(:,i)~=0);
    subplot(2*nTemplates, 1, nTemplates+i)
    hold on
    lineColors = hsv(nTraces);
    if (sum(failuresBin) || isfield(artifact, 'estimate')) && nSpikes
        for k = 1:nTraces %each trace
            if latenciesFull(k,i)
                for m = 1:nChannels % each electrode
                    if m == 1 || sum(failuresBin)
                        offset = (m-1)*28;
                        current = plot(1+offset:26+offset, squeeze(ei(k,m,:)));
                        set(findobj(current,'Type','line'),'Color',lineColors(k,:))
                    end
                end
            end
        end
    end
    for j = 1:nChannels
        offset = (j-1)*28;
        plot(1+offset:26+offset, squeeze(targetEi(j,targetEiMinPos-10:targetEiMinPos+15)),'k')
    end
    hold off
    set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
    
    
    ylimits = get(gca, 'ylim');
    for j = 1:nChannels
        text((28*(j-1)+13), (ylimits(2)-ylimits(1))*0.95+ylimits(1), num2str(channelsToUse(j)))
    end

    title(['neuron ' num2str(neuronIDs(i)) ' (', num2str(nSpikes),')'])
    
end

set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize=26;
xLeft = (22-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])


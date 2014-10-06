function figHandle = summaryPlotGeneratorBrief(dataTraces, centerChannel, channelRadius, data,...
    artifact, template, results, shift)

% a modified version of summaryPlotGenerator that only plots the basics (raw data, and data after
% subtraction of mean artifact trace)
%
% used to generate plot for ARVO 2009 poster
%



dataPath = data.path;
patternNumber = data.patternNumber;
movieNumber = data.movieNumber;

actualArtifact = artifact.actual;
actualArtifactInit = artifact.actualInit;
artifactMean = artifact.mean;
modelType = artifact.modelType;
if isfield(artifact, 'estimate')
    estArtifact = artifact.estimate;
end

templates = template.vectors;
pathToEi = template.pathToEi;
neuronIDs = template.neuronIDs;

latenciesInit = results.latenciesInit;
latenciesFull = results.latenciesFull;
successesInit = results.successesInit;
successesFull = results.successesFull;


nTraces  = size(dataTraces, 1);
nSamples = size(dataTraces, 3);
nTemplates = length(templates);


ymin1 = min(min(squeeze(dataTraces(:,centerChannel,:))));
ymax1 = max(max(squeeze(dataTraces(:,centerChannel,:))));
ymin1 = -500;
ymax1 = 50*ceil(ymax1/50);

ymin2 = -100;
ymax2 = 100;

xPlotRange = 50;

figure('position', [100 100 400 500])
subplot(2,1,1)
hold on

failuresBin = zeros(1, nTraces);
nSuccesses = 0;
if sum(reshape(successesFull,[],1))
    successesToPlotFull = [];
    failuresToPlotFull = [];
    for i = 1:nTraces
        if any(successesFull(i,:))
            successesToPlotFull = [successesToPlotFull squeeze(dataTraces(i,centerChannel,:))]; %#ok<AGROW>
            nSuccesses = nSuccesses + 1;
        else
            failuresToPlotFull = [failuresToPlotFull squeeze(dataTraces(i,centerChannel,:))]; %#ok<AGROW>
            failuresBin(i) = 1;
        end
    end
    if ~isempty(failuresToPlotFull)
        current = plot((1:nSamples)/20, failuresToPlotFull);
        set(findobj(current,'Type','line'),'Color',[0 0 0])
    end
    plot((1:nSamples)/20, successesToPlotFull, 'r-')
end
hold off

title('raw data', 'fontsize', 20)
xlabel('time (ms)', 'fontsize', 20)
set(gca, 'ytick', [])
set(gca, 'xlim', [1/20 xPlotRange/20], 'ylim', [ymin1 ymax1], 'fontsize', 18)


%%%%%%%%%%%%%% plots traces with spikes after artifact subtraction
subplot(2,1,2)

hold on
%nSuccesses = 0;
for i = 1:nTraces
    if max(successesInit(i,:))
        %nSuccesses = nSuccesses + 1;
        plot((1:nSamples)/20, squeeze(dataTraces(i, centerChannel, :)) - actualArtifactInit,'r')
    end
end
%    templateTimeVector = (medLatencyInit - templateMinPos + 1):nSamples;
%    plot(templateTimeVector, template(1:length(templateTimeVector)), 'k')
hold off
title('after subtraction of mean artifact trace', 'fontsize', 20)
set(gca, 'xlim', [1 xPlotRange]/20, 'ylim', [ymin2 ymax2], 'fontsize', 18)
xlabel('time (ms)', 'fontsize', 20)
set(gca, 'ytick', [])


keyboard


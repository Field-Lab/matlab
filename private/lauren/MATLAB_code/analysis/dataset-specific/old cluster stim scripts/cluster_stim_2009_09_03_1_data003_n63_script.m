clear all

plotLogScale = 1;
constrainSlopes = 0;

pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data003';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 63;

pElec = 7;

patternNos = 1:253;
movieNosRange = 1:373;

%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)


%%
[thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID,...
    relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 0,...
    'constrainSlopes', constrainSlopes, 'threshPlotLims', [0.07 0.15], 'curvePlotXLims', [0.05 1],...
    'redoConstrainedFitting', 0, 'excludePatterns', [], 'recalcAll', 0, 'insetYLim', [0 0.3]);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)
keyboard

generateSecondaryCombPlots(details.secondaryCombsIncluded, thresholds, threshStds, estPrimThresh, details.actualRelAmps)


%% plotting model and data for electrode pair

figure('position', [100 100 400 300]*2)
axes('units', 'pixels', 'position', [20 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 1,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 .3])

axes('units', 'pixels', 'position', [20 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 2,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 0.3])

axes('units', 'pixels', 'position', [20 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 3,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 0.3])

axes('units', 'pixels', 'position', [220 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 4,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 0.3])

axes('units', 'pixels', 'position', [220 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 5,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 0.3])

axes('units', 'pixels', 'position', [220 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 6,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-.3 .3], 'ylim', [0 0.3])


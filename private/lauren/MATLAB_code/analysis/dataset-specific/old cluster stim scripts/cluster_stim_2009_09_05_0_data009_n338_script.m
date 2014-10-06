clear all

plotLogScale = 1;
constrainSlopes = 0;
plotInconsistentFits = true;

pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-05-0/data009';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 338;

pElec = 28;

patternNos = 266:530;
movieNosRange = 1:221;

excludePatterns = [297 298]; %falls outside limits set by secondary-alone thresholds (2x anodal stimulation)

%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)

%%
[thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID,...
    relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 0,...
    'constrainSlopes', constrainSlopes, 'threshPlotLims', [0.07 0.15], 'curvePlotXLims', [0.14 1],...
    'redoConstrainedFitting', 0, 'recalcAll', 0, 'plotInconsistentFits', plotInconsistentFits, 'excludePatterns', excludePatterns);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

keyboard

generateSecondaryCombPlots(details.secondaryCombsIncluded, thresholds, threshStds, estPrimThresh, details.actualRelAmps)

%%

showSecAloneThresh = true;
fitWithinSecAloneThresh = true;

figure('position', [100 100 400 300]*2)
axes('units', 'pixels', 'position', [20 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 1,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])

axes('units', 'pixels', 'position', [20 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 2,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])

axes('units', 'pixels', 'position', [20 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 3,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-2 1.2], 'ylim', [0 1.6])

axes('units', 'pixels', 'position', [220 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 4,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])

axes('units', 'pixels', 'position', [220 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 5,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])

axes('units', 'pixels', 'position', [220 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 6,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])
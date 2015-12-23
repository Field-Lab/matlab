clear all

plotLogScale = 1;
constrainSlopes = 0;

pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data014_015';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 91;

pElec = 7;

patternNos = 1:253;
movieNosRange = 1:320;

%%
[thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID,...
    relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 1,...
    'constrainSlopes', constrainSlopes, 'threshPlotLims', [0.2 0.4], 'curvePlotXLims', [0.15 0.6],...
    'redoConstrainedFitting', 0, 'excludePatterns', [], 'recalcAll', 0);

keyboard


%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

generateSecondaryCombPlots(details.secondaryCombsIncluded, thresholds, threshStds, estPrimThresh, details.actualRelAmps)


%% plot model vs. data

figure('position', [100 100 400 300]*2)
axes('units', 'pixels', 'position', [20 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 1)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])

axes('units', 'pixels', 'position', [20 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 2)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])

axes('units', 'pixels', 'position', [20 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 3)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])

axes('units', 'pixels', 'position', [220 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 4)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])

axes('units', 'pixels', 'position', [220 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 5)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])

axes('units', 'pixels', 'position', [220 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 6)
set(gca, 'xlim', [-.8 .8], 'ylim', [0 .8])
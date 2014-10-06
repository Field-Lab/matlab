clear all

plotLogScale = 1;
constrainSlopes = 0;

showSecAloneThresh = true;
fitWithinSecAloneThresh = true;



pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data014_015';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 409; 

pElec = 28;

patternNos = 266:530;
movieNosRange = 1:320;

%%

% [thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID,...
%     relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 1,...
%     'constrainSlopes', constrainSlopes, 'threshPlotLims', [0.1 1], 'curvePlotXLims', [0.15 1],...
%     'redoConstrainedFitting', 0, 'excludePatterns', [], 'recalcAll', 0, 'insetYLim', [0 1.2]);

[thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID,...
    relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 1,...
    'constrainSlopes', constrainSlopes, 'threshPlotLims', [0.1 1], 'curvePlotXLims', [0.15 1],...
    'redoConstrainedFitting', 0, 'excludePatterns', [], 'recalcAll', 0, 'insetYLim', [0 1.2], 'excludePatterns', [286 292 298]);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

generateSecondaryCombPlots(details.secondaryCombsIncluded, thresholds, threshStds, estPrimThresh, details.actualRelAmps)

%% plot model vs. data

figure('position', [100 100 400 320]*2)
axes('units', 'pixels', 'position', [20 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 1,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(1))])

axes('units', 'pixels', 'position', [20 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 2,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(2))])

axes('units', 'pixels', 'position', [20 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 3,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(3))])

axes('units', 'pixels', 'position', [220 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 4,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(4))])

axes('units', 'pixels', 'position', [220 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 5,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(5))])

axes('units', 'pixels', 'position', [220 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 6,...
    'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
set(gca, 'xlim', [-1.5 1.5], 'ylim', [0 1.5])
title(['secondary electrode ' num2str(details.sElecs(6))])
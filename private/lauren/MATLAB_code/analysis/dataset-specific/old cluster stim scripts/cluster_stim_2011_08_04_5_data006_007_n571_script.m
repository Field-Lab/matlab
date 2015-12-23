clear all

plotLogScale = 1; %only affects plots, not curve fits
constrainSlopes = 0;
showSecAloneThresh = true;
fitWithinSecAloneThresh = true;
plotInconsistentFits = false;

pathToData = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data006_007';

relAmps = [-1.5 -1 -0.5 -0.25 0.25 0.5 1 1.5];
neuronID = 571;

pElec = 43;
patternNos = [1:72 290 292:303];
excludePatterns = [];


%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID,...
    'xPlotLim', [-2 2], 'yPlotLim', [1 3], 'checkForPAloneArt', true,...
    'fitWithinSecAloneThresh', fitWithinSecAloneThresh, 'showSecAloneThresh', showSecAloneThresh)

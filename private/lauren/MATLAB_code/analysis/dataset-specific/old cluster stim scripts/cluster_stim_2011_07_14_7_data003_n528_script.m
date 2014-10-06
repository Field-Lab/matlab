clear all

plotLogScale = 1; %only affects plots, not curve fits
constrainSlopes = 0;
showSecAloneThresh = true;
fitWithinSecAloneThresh = true;
plotInconsistentFits = false;

pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 528;

pElec = 41;
patternNos = [304:339 593:2:605];
excludePatterns = [];


%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID,...
    'xPlotLim', [-1 1], 'yPlotLim', [0.3 1.3], 'checkForPAloneArt', true)

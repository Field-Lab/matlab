clear all

plotLogScale = 1; %only affects plots, not curve fits
constrainSlopes = 0;
showSecAloneThresh = true;
fitWithinSecAloneThresh = true;
plotInconsistentFits = false;

pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 632;

pElec = 43;
patternNos = [607:642 896:2:908];
excludePatterns = [];


%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID, 'xPlotLim', [-0.5 0.5], 'yPlotLim', [0.5 1], 'checkForPAloneArt', true)

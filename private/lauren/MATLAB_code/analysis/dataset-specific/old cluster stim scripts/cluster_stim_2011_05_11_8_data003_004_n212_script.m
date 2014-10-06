clear all

plotLogScale = 1; %only affects plots, not curve fits
constrainSlopes = 0;
showSecAloneThresh = true;
fitWithinSecAloneThresh = true;
plotInconsistentFits = false;

pathToData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data003_004';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 212;

pElec = 14;
patternNos = [1:36 290:2:302];
excludePatterns = [];
%movieNosRange = 1:126;


%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)

%%
[thresholds threshStds details] = generateBasicClusterStimPlotPairsOnly(pathToData, patternNos,...
    pElec, neuronID, 'xPlotLim', [-0.4 0.4], 'yPlotLim', [0.4 0.8], 'checkForPAloneArt', true);



% [thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec,...
%     neuronID, relAmps, 'plotLogScale', plotLogScale, 'plotEstimatedPrimThresh', 0, 'constrainSlopes', constrainSlopes,...
%     'threshPlotLims', [0.25 0.8], 'curvePlotXLims', [0.15 1], 'redoConstrainedFitting', 0,...
%     'recalcAll', 0, 'plotInconsistentFits', plotInconsistentFits, 'excludePatterns', excludePatterns);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

%generateSecondaryCombPlots(details.secondaryCombsIncluded, thresholds, threshStds, estPrimThresh, details.actualRelAmps)

% 
% %% plotting model and data for electrode pair
% 
% slopes = zeros(6,1);
% 
% figure('position', [100 100 400 300]*2)
% for i = 1:6
%     if i<4 %first column
%         axes('units', 'pixels', 'position', [20 320-100*i 160 80]*2)
%     else
%         axes('units', 'pixels', 'position', [220 320-100*(i-3) 160 80]*2)
%     end
%     slopes(i) = plotLinearityDataVsModel(thresholds, threshStds, details, i,...
%         'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh);
%     set(gca, 'xlim', [-1 1], 'ylim', [0 1])
% end
% 
% 
% %ylabel('primary electrode current')
% %xlabel('secondary electrode current')
% 
% %% plotting "electrical receptive field"
% markerScale = 100;
% plotElectricRF(details.sElecs, slopes, markerScale, pElec)
% 
% 
% %%
% 
% showSecAloneThresh = true;
% fitWithinSecAloneThresh = true;
% 
% figure('position', [100 100 400 300]*2)
% axes('units', 'pixels', 'position', [20 220 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 1,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])
% 
% axes('units', 'pixels', 'position', [20 120 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 2,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])
% 
% axes('units', 'pixels', 'position', [20 20 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 3,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-2 1.2], 'ylim', [0 1.6])
% 
% axes('units', 'pixels', 'position', [220 220 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 4,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])
% 
% axes('units', 'pixels', 'position', [220 120 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 5,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])
% 
% axes('units', 'pixels', 'position', [220 20 160 80]*2)
% plotLinearityDataVsModel(thresholds, threshStds, details, 6,...
%     'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh)
% set(gca, 'xlim', [-1.25 1.25], 'ylim', [0 1.25])






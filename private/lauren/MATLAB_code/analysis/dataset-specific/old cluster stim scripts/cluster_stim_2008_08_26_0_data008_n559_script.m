plotLogScale = 1;
constrainSlopes = 0;

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008';

relAmps = [-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1];
neuronID = 559; 

pElec = 44;
%patternNos = [221:318];
patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317]; %only patterns with negative primary
movieNosRange = 1:126;

%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)

%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID, 'xPlotLim', [-1 1], 'yPlotLim', [0 1],...
    'checkForPAloneArt', true)

[thresholds, threshStds, estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec,...
    neuronID, relAmps, 'plotLogScale', plotLogScale, 'constrainSlopes', constrainSlopes,...
    'threshPlotLims', [0.40 0.65], 'curvePlotXLims', [0.25 1.2], 'redoConstrainedFitting', 0, 'recalcAll', 0);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

%% plotting model and data for electrode pairs

figure('position', [100 100 400 300]*2)
axes('units', 'pixels', 'position', [20 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 1)
set(gca, 'xlim', [-1 1], 'ylim', [0 1])

axes('units', 'pixels', 'position', [20 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 2)
set(gca, 'xlim', [-1 1], 'ylim', [0 1])

axes('units', 'pixels', 'position', [20 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 3)
set(gca, 'xlim', [-1 1], 'ylim', [0 1])

axes('units', 'pixels', 'position', [220 220 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 4)
set(gca, 'xlim', [-1 1], 'ylim', [0 1])

axes('units', 'pixels', 'position', [220 120 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 5)
set(gca, 'xlim', [-1 1], 'ylim', [0 1])

axes('units', 'pixels', 'position', [220 20 160 80]*2)
plotLinearityDataVsModel(thresholds, threshStds, details, 6)
set(gca, 'xlim', [-1 1], 'ylim', [0 1.5])


%xlabel('primary electrode current')
%ylabel('secondary electrode current')

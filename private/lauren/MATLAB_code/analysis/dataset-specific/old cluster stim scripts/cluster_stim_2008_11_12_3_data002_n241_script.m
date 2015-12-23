clear all

plotLogScale = 1;
constrainSlopes = 0;

showSecAloneThresh = true;
fitWithinSecAloneThresh = true;

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002';

relAmps = [-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1];
neuronID = 241; 

pElec = 21;
%patternNos = [221:318];
patternNos = [1:6 13:18 25:30 37:42 49:54 61:66 73:78 85:90 97:110];
movieNosRange = 1:126;

%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)

%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID, 'xPlotLim', [-1.5 1.5], 'yPlotLim', [1.25 2.75],...
    'checkForPAloneArt', true)

[thresholds, threshStds, estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec,...
    neuronID, relAmps, 'plotLogScale', plotLogScale, 'constrainSlopes', constrainSlopes,...
    'threshPlotLims', [0.40 0.65], 'curvePlotXLims', [1 4], 'redoConstrainedFitting', 0,...
    'recalcAll', 0, 'erfStartParams', [1 -2.5], 'insetYLim', [0 4]);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)

%% plotting model and data for electrode pairs

slopes = zeros(6,1);
figure('position', [100 100 400 300]*2)
for i = 1:6
    if i<4 %first column
        axes('units', 'pixels', 'position', [20 320-100*i 160 80]*2)
    else
        axes('units', 'pixels', 'position', [220 320-100*(i-3) 160 80]*2)
    end
    slopes(i) = plotLinearityDataVsModel(thresholds, threshStds, details, i,...
        'showSecAloneThresh', showSecAloneThresh, 'fitWithinSecAloneThresh', fitWithinSecAloneThresh);
    set(gca, 'xlim', [-4.5 4.5], 'ylim', [0 4.5])
end


%% plotting "electrical receptive field"
markerScale = 100;
plotElectricRF(details.sElecs, slopes, markerScale, pElec)

















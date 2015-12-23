plotLogScale = 1;
constrainSlopes = 0;
showSecAloneThresh = true;
fitWithinSecAloneThresh = true;

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';

relAmps = [-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1];
neuronID = 800; 

pElec = 49;
%patternNos = [221:318];
patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317:330]; %only patterns with negative primary
movieNosRange = 1:126;

%bootstrapNDrawsCompare(pathToData, patternNos, neuronID, 'recalcAll', 0, 'bootstrapReps', 100)

%%
generateBasicClusterStimPlotPairsOnly(pathToData, patternNos, pElec, neuronID, 'xPlotLim', [-0.4 0.4], 'yPlotLim', [0.3 0.7],...
    'checkForPAloneArt', true)

[thresholds threshStds estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec,...
    neuronID, relAmps, 'plotLogScale', plotLogScale, 'constrainSlopes', constrainSlopes,...
    'threshPlotLims', [0.25 0.8], 'curvePlotXLims', [0.15 1.5], 'redoConstrainedFitting', 0,...
    'recalcAll', 0);

%checkSystematicResids(pathToData, patternNos, neuronID, movieNosRange)


%% plotting model and data for electrode pair

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
    set(gca, 'xlim', [-1 1], 'ylim', [0 1])
end


%ylabel('primary electrode current')
%xlabel('secondary electrode current')

%% plotting "electrical receptive field"
markerScale = 100;
plotElectricRF(details.sElecs, slopes, markerScale, pElec)

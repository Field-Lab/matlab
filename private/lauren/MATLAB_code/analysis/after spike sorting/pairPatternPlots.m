function pairPatternPlots(thresholds, threshStds, details, varargin)
%function [slopes yInts excludedFromFits] = pairPatternPlots(thresholds, threshStds, details, varargin)


p = inputParser;

p.addRequired('thresholds', @isnumeric)
p.addRequired('threshStds', @isnumeric)
p.addRequired('details', @isstruct)

p.addParamValue('xPlotLim', [0 1.5], @isnumeric)
p.addParamValue('yPlotLim', [0 1], @isnumeric)
p.addParamValue('checkForPAloneArt', false, @islogical)
p.addParamValue('plotAxes',[], @iscell) %cell array of axes handles for each secondary electrode plot (in order given by getCluster)
p.addParamValue('showSecAloneThresh', false, @islogical)
p.addParamValue('plotUnfitPairs', true, @islogical) %whether to plot data points that aren't used in fitting (in any region)

p.parse(thresholds, threshStds, details, varargin{:})

checkForPAloneArt = p.Results.checkForPAloneArt;
plotAxes = p.Results.plotAxes;
xPlotLim = p.Results.xPlotLim;
yPlotLim = p.Results.yPlotLim;


%% separates primary-alone data into 6 subsets and calculates
% threshold based on subset of data, to see if there is any systematic
% deviation of p-alone from linearity

%NEEDS TO BE UPDATED TO WORK WITH CURRENT CODE

% if checkForPAloneArt
%     params.xPlotLim = xPlotLim;
%     params.yPlotLim = yPlotLim;
%     params.showSecAloneThresh = p.Results.showSecAloneThresh;
%     plotCheckForPAloneArt(pAloneBin, details.pAmpsAll, neuronID, patternNos, thresholds, threshStds, details, params);
% end


%% plots of thresholds and comparison to linear model

if isempty(plotAxes)
    figure('position', [100 100 400 300]*2)
    axesPos = {[20 220 160 80]*2, [20 120 160 80]*2, [20 20 160 80]*2, [220 220 160 80]*2, [220 120 160 80]*2, [220 20 160 80]*2};
end


for jj = 1:6
    if isempty(plotAxes)
        axes('units', 'pixels', 'position', axesPos{jj})
    else
        axes(plotAxes{jj})
    end

    plotLinearityDataVsModel(thresholds, threshStds, details, jj,...
        'showSecAloneThresh', p.Results.showSecAloneThresh,...
        'xPlotLim', xPlotLim, 'yPlotLim', yPlotLim, 'plotUnfitPairs', p.Results.plotUnfitPairs, 'shadeRegType', 'bootstrapSD');
    
    set(gca, 'xlim', xPlotLim, 'ylim', yPlotLim, 'box', 'on')
    title(['electrode ' num2str(details.sElecs(jj))])
    xlabel('secondary current amplitude')
    ylabel('primary current amplitude')
    
end


end


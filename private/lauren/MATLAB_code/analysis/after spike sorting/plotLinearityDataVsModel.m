function plotLinearityDataVsModel(thresholds, threshStds, details, sInd, varargin)

%details: output of extractPatternStimAnalysis.m

p = inputParser;

p.addRequired('thresholds', @isnumeric)
p.addRequired('threshStds', @isnumeric)
p.addRequired('details', @isstruct)
p.addRequired('sInd', @isnumeric)

p.addParamValue('showSecAloneThresh', false, @islogical)
p.addParamValue('showSecAloneThreshMin', false, @islogical)
p.addParamValue('xPlotLim', [], @isnumeric)
p.addParamValue('yPlotLim', [], @isnumeric)
p.addParamValue('plotUnfitPairs', true, @islogical) %whether to plot data points that aren't used in fitting (in any region)
p.addParamValue('shadeRegType', 'curveSD', @ischar) %"curveSD" or "bootstrapSD" -- determines whether gray region represents
% 1) dynamic range of response curves or 2) standard deviation of threshold
% estimates from bootstrap analysis

%p.addParamValue('sElecConnected', [], @iscell) %this is now in the 'details' struct


p.parse(thresholds, threshStds, details, sInd, varargin{:})

showSecAloneThresh = p.Results.showSecAloneThresh;
showSecAloneThreshMin = p.Results.showSecAloneThreshMin;
xPlotLim = p.Results.xPlotLim;
yPlotLim = p.Results.yPlotLim;
plotUnfit = p.Results.plotUnfitPairs;

if strcmpi(p.Results.shadeRegType, 'curveSD')
    shadeBS = false;
elseif strcmpi(p.Results.shadeRegType, 'bootstrapSD')
    shadeBS = true;
else
    error('unrecognized parameter value for ''shadeRegType''')
end

%sConnected = p.Results.sElecConnected;


nPatterns = length(thresholds);

details.pAloneBin = cast(details.pAloneBin, 'logical');

if isfield(details, 'sConnected')
    sConnected = details.sConnected;
else
    sConnected = [];
end


%% extracts data

hold on

xToPlot = [];
yToPlot = [];
excludedFromFits = logical([]);
xToPlotBSSDLow = []; %bootstrap standard deviations
yToPlotBSSDLow = [];
xToPlotBSSDHigh = [];
yToPlotBSSDHigh = [];
xToPlotSDLow = []; %standard deviation of cumulative gaussian fit
yToPlotSDLow = [];
xToPlotSDHigh = [];
yToPlotSDHigh = [];
relAmps = [];
yToPlotConnected = [];
toPlotPInd = [];

for i = 1:nPatterns
    % primary alone or electrode pairs
    if (any(details.sAmps{i}(sInd,:)) && sum(any(details.sAmps{i},2) ~= 0) == 1 && any(details.pAmps{i})) &&...
            ~isinf(thresholds(i)) || (details.pAloneBin(i) && details.pAmps{i}(1) < 0)
        w = details.curveWidths(i);
        
        relAmp = mean(details.actualRelAmps{i}(sInd,:)); %since it can vary slightly across amplitudes, just take the mean
        
        xToPlot = [xToPlot thresholds(i)*relAmp]; %#ok<AGROW>
        yToPlot = [yToPlot thresholds(i)]; %#ok<AGROW>
        
        %details.unfitPairs signifies pairs that aren't used in the
        %fitting in any region because the threshold is infinite, the
        %secondary currents sometimes round to zero, or the pattern is
        %specifically marked as bad in 'cell_list_cluster_stim'
        excludedFromFits = [excludedFromFits details.unfitPairs(i)];
        toPlotPInd = [toPlotPInd i];
        
        xToPlotBSSDLow = [xToPlotBSSDLow (thresholds(i)-threshStds(i))*relAmp];
        yToPlotBSSDLow = [yToPlotBSSDLow (thresholds(i)-threshStds(i))];
        xToPlotBSSDHigh = [xToPlotBSSDHigh (thresholds(i)+threshStds(i))*relAmp];
        yToPlotBSSDHigh = [yToPlotBSSDHigh (thresholds(i)+threshStds(i))];
        
        xToPlotSDLow = [xToPlotSDLow thresholds(i)*relAmp*(1-1/w)]; %#ok<AGROW>
        yToPlotSDLow = [yToPlotSDLow thresholds(i)*(1-1/w)];     %#ok<AGROW>
        xToPlotSDHigh = [xToPlotSDHigh thresholds(i)*relAmp*(1+1/w)]; %#ok<AGROW>
        yToPlotSDHigh = [yToPlotSDHigh thresholds(i)*(1+1/w)];  %#ok<AGROW>
        
        relAmps = [relAmps relAmp]; %#ok<AGROW>
        
        if details.pAloneBin(i)
            pAloneThresh = thresholds(i);
        end

    % secondary electrode alone
    elseif showSecAloneThresh && details.sAloneBin(i)
        if details.sAmps{i}(sInd,1) < 0 %corresponds to secondary-alone cathodal
            sAlone(1).thresh = thresholds(i);
            if ~isinf(thresholds(i))
                w = details.curveWidths(i);
                               
                sAlone(1).xBSSDLow =  thresholds(i)-threshStds(i);
                sAlone(1).yBSSDLow =  0;
                sAlone(1).xBSSDHigh = thresholds(i)+threshStds(i);
                sAlone(1).yBSSDHigh = 0;
                
                sAlone(1).xSDLow =  thresholds(i)*(1-1/w);
                sAlone(1).ySDLow =  0;
                sAlone(1).xSDHigh = thresholds(i)*(1+1/w);
                sAlone(1).ySDHigh = 0;
            end

            if isfield(details, 'threshMins')
                sAlone(1).threshMins(1) = details.threshMins(i);
            else
                sAlone(1).threshMins(1) = 0;
            end
        elseif details.sAmps{i}(sInd,1) > 0 %corresponds to secondary-alone anodal
            sAlone(2).thresh = -thresholds(i);
            if ~isinf(thresholds(i))
                w = details.curveWidths(i);
                
                sAlone(2).xBSSDLow =  -1*(thresholds(i)-threshStds(i));
                sAlone(2).yBSSDLow =  0;
                sAlone(2).xBSSDHigh = -1*(thresholds(i)+threshStds(i));
                sAlone(2).yBSSDHigh = 0;
                
                sAlone(2).xSDLow =  -1*thresholds(i)*(1-1/w);
                sAlone(2).ySDLow =  0;
                sAlone(2).xSDHigh = -1*thresholds(i)*(1+1/w);
                sAlone(2).ySDHigh = 0;
            end
            if isfield(details, 'threshMins')
                sAlone(2).threshMins(2) = -details.threshMins(i);
            else
                sAlone(2).threshMins(2) = 0;
            end
        end
    end
    
    
    % plot thresholds for connected but zero secondary current secondary-primary pairs
    if ~isempty(sConnected)
        if sConnected{i}(sInd) && ~any(details.sAmps{i}(sInd,:))
            yToPlotConnected = [yToPlotConnected thresholds(i)];
        end
    end
end


%% plotting

regCol = lines(3);
regColShade = 1 - (1-regCol)*0.3;


mBound = details.mValidBound(:,sInd);

%determine which data points lie in which regions
% region 1 = region with primary-alone threshold
regBin{2} = false(size(xToPlot)); %secondary positive (cathodal) region
regBin{3} = false(size(xToPlot)); %secondary negative (anodal) region
if mBound(1)
    regBin{2}(xToPlot>0 & abs(yToPlot./xToPlot) < mBound(1)) = true;
end
if mBound(2)
    regBin{3}(xToPlot<0 & abs(yToPlot./xToPlot) < abs(mBound(2))) = true;
end
regBin{1} = ~(regBin{2} | regBin{3});


hold on

if plotUnfit
    excludeFromReg = false(size(xToPlot));
else
    excludeFromReg = excludedFromFits;
end

%plot shaded regions (either response curve width or
% threshold standard deviation from bootstrap analysis)
for ii = 1:3
    if shadeBS % shading corresponds to thresholds SD from bootstrap analysis
        xShadeRegLow{ii} =  xToPlotBSSDLow(regBin{ii}  & ~excludeFromReg); %#ok<*AGROW>
        xShadeRegHigh{ii} = xToPlotBSSDHigh(regBin{ii} & ~excludeFromReg);
        yShadeRegLow{ii} =  yToPlotBSSDLow(regBin{ii}  & ~excludeFromReg);
        yShadeRegHigh{ii} = yToPlotBSSDHigh(regBin{ii} & ~excludeFromReg);
        
        relAmpsReg{ii} = relAmps(regBin{ii} & ~excludeFromReg);
        
        %include s-alone in region
        if ii == 2 && ~isempty(xShadeRegLow{ii});
            xShadeRegLow{ii} =  [xShadeRegLow{ii}  sAlone(1).xBSSDLow];
            xShadeRegHigh{ii} = [xShadeRegHigh{ii} sAlone(1).xBSSDHigh];
            yShadeRegLow{ii} =  [yShadeRegLow{ii}  sAlone(1).yBSSDLow];
            yShadeRegHigh{ii} = [yShadeRegHigh{ii} sAlone(1).yBSSDHigh];
            relAmpsReg{ii} =    [relAmpsReg{ii}    inf];
        elseif ii == 3 && ~isempty(xShadeRegLow{ii});
            xShadeRegLow{ii} =  [xShadeRegLow{ii}  sAlone(2).xBSSDLow];
            xShadeRegHigh{ii} = [xShadeRegHigh{ii} sAlone(2).xBSSDHigh];
            yShadeRegLow{ii} =  [yShadeRegLow{ii}  sAlone(2).yBSSDLow];
            yShadeRegHigh{ii} = [yShadeRegHigh{ii} sAlone(2).yBSSDHigh];
            relAmpsReg{ii} =    [relAmpsReg{ii}    -inf];
        end
    else % shading corresponds SD of cumulative Gaussian fit to response curve (dynamic range)
        xShadeRegLow{ii} =  xToPlotSDLow(regBin{ii}  & ~excludeFromReg);
        xShadeRegHigh{ii} = xToPlotSDHigh(regBin{ii} & ~excludeFromReg);
        yShadeRegLow{ii} =  yToPlotSDLow(regBin{ii}  & ~excludeFromReg);
        yShadeRegHigh{ii} = yToPlotSDHigh(regBin{ii} & ~excludeFromReg);
        
        relAmpsReg{ii} = relAmps(regBin{ii} & ~excludeFromReg);
        
        %include s-alone in region
        if ii == 2 && ~isempty(xShadeRegLow{ii});
            xShadeRegLow{ii} =  [xShadeRegLow{ii}  sAlone(1).xSDLow];
            xShadeRegHigh{ii} = [xShadeRegHigh{ii} sAlone(1).xSDHigh];
            yShadeRegLow{ii} =  [yShadeRegLow{ii}  sAlone(1).ySDLow];
            yShadeRegHigh{ii} = [yShadeRegHigh{ii} sAlone(1).ySDHigh];
            relAmpsReg{ii} =    [relAmpsReg{ii}    inf];
        elseif ii == 3 && ~isempty(xShadeRegLow{ii});
            xShadeRegLow{ii} =  [xShadeRegLow{ii}  sAlone(2).xSDLow];
            xShadeRegHigh{ii} = [xShadeRegHigh{ii} sAlone(2).xSDHigh];
            yShadeRegLow{ii} =  [yShadeRegLow{ii}  sAlone(2).ySDLow];
            yShadeRegHigh{ii} = [yShadeRegHigh{ii} sAlone(2).ySDHigh];
            relAmpsReg{ii} =    [relAmpsReg{ii}    -inf];
        end
        
    end
    
    [~, sortInd] = sort(relAmpsReg{ii});
    [~, sortIndRev] = sort(relAmpsReg{ii}, 2, 'descend');
    
    fillPolyX = [xShadeRegLow{ii}(sortInd) xShadeRegHigh{ii}(sortIndRev)];
    fillPolyY = [yShadeRegLow{ii}(sortInd) yShadeRegHigh{ii}(sortIndRev)];

    fill(fillPolyX, fillPolyY, regColShade(ii,:), 'edgecolor', 'none')
    
    plot(xToPlot(regBin{ii} & ~excludedFromFits), yToPlot(regBin{ii} & ~excludedFromFits), 'o', 'MarkerEdgeColor', regCol(ii,:))
    
    if plotUnfit
        plot(xToPlot(regBin{ii} & excludedFromFits), yToPlot(regBin{ii} & excludedFromFits), 'x', 'MarkerEdgeColor', regCol(ii,:));
        
        %plot thresholds for connected but 0-secondary-current combinations
        if ~isempty(yToPlotConnected)
            plot(zeros(length(yToPlotConnected)), yToPlotConnected, 'k*')
        end
    end
    
    %plots bootstrap standard deviations if they aren't represented by
    %shaded region
    if ~shadeBS
        for jj = 1:length(xToPlot)
            if ~excludedFromFits(jj) && regBin{ii}(jj)
                plot([xToPlotBSSDLow(jj) xToPlotBSSDHigh(jj)], [yToPlotBSSDLow(jj) yToPlotBSSDHigh(jj)], '-', 'color', regCol(ii,:))
            end
        end
    end
end


%plot secondary-alone thresholds
if showSecAloneThresh
    if ~isinf(sAlone(1).thresh) %cathodal (positive x axis)
        plot(sAlone(1).thresh, 0, 'o', 'MarkerEdgeColor', regCol(2,:))
        plot([sAlone(1).xBSSDLow sAlone(1).xBSSDHigh], [sAlone(1).yBSSDLow sAlone(1).yBSSDHigh], '-', 'color', regCol(2,:))
    end
    if ~isinf(sAlone(2).thresh) %anodal
        plot(sAlone(2).thresh, 0, 'o', 'MarkerEdgeColor', regCol(3,:))
        plot([sAlone(2).xBSSDLow sAlone(2).xBSSDHigh], [sAlone(2).yBSSDLow sAlone(2).yBSSDHigh], '-', 'color', regCol(3,:))
    end
    if showSecAloneThreshMin
        if sAloneThreshMins(1) %zero signifies none
            plot(sAlone(1).threshMins, 0, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
        end
        if sAloneThreshMins(2) %zero signifies none
            plot(sAlone(2).ThreshMins, 0, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
        end
    end
end


%%
model = details.multiRegModel(sInd);

%%% plots model line(s) %%%
if mBound(1)
    alpha = [0 model.posInt(2)];
    plot(model.tsPos - alpha*model.sLamPos, alpha, '-', 'color', regCol(2,:))
end
if mBound(2)
    alpha = [0 model.negInt(2)];
    plot(-model.tsNeg - alpha*model.sLamNeg, alpha, '-', 'color', regCol(3,:))
end

if mBound(1) && mBound(2)
    alpha = [model.negInt(1)         model.posInt(1)];
elseif mBound(1)
    alpha = [min(xToPlot(regBin{1})) model.posInt(1)];
elseif mBound(2)
    alpha = [model.negInt(1)         max(xToPlot(regBin{1}))];
else
    alpha = [min(xToPlot(regBin{1})) max(xToPlot(regBin{1}))];
end
plot(alpha, model.tp - alpha*model.pLam, '-', 'color', regCol(1,:))


%plot boundaries between regions
if mBound(1)
    plot([0 1], [0 mBound(1)], 'k--')
end
if mBound(2)
    plot([0 -1], [0 -mBound(2)], 'k--')
end

hold off
if ~isempty(xPlotLim) && ~isempty(yPlotLim)
    set(gca, 'xlim', xPlotLim, 'ylim', yPlotLim)
end


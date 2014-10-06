function [r_sq unfitTripBin badCompPairBin] = generateSecondaryCombPlots_new(thresholds, threshStds, details, varargin)

%
%
%

p = inputParser;

p.addRequired('thresholds', @isnumeric)
p.addRequired('threshStds', @isnumeric)
p.addRequired('details', @isstruct)

%p.addParamValue('xPlotLim', [0 1.5], @isnumeric)
%p.addParamValue('yPlotLim', [0 1], @isnumeric)
p.addParamValue('estPrimThresh', [], @isnumeric)% estimated primary-alone threshold to use in place of actual measured threshold
p.addParamValue('obsVExpAxes', [], @ishandle)
p.addParamValue('plotUnfitted', true, @islogical) %whether to still plot the points that aren't included in the linear fit
p.addParamValue('obsVExpColor', [], @isnumeric) %if specified, all points in observed vs expected plot are plotted in this color
p.addParamValue('useFitPAloneThresh', false, @islogical) %if true, uses mean of y-intercepts from electrode pair linear fits as primary-alone threshold, rather than actual measured primary-alone threshold

% whether to use actual component pair thresholds shifts for predicted triplet shift (false)
%or to predict triplet shifts from slopes of best-fit line from component
%pairs (true)
p.addParamValue('fullModel', false, @islogical)

p.parse(thresholds, threshStds, details, varargin{:})

%xPlotLim = p.Results.xPlotLim;
%yPlotLim = p.Results.yPlotLim;
estPrimThresh = p.Results.estPrimThresh;
fullModel = p.Results.fullModel;
obsVExpAxes = p.Results.obsVExpAxes;
plotUnfitted = p.Results.plotUnfitted;
obsVExpColor = p.Results.obsVExpColor;
useFitPAloneThresh = p.Results.useFitPAloneThresh;

%% sort out which component pairs are associated with each triplet
nPatterns = length(thresholds);


%first, determine which patterns correspond to electrode pairs
pairFlags = false(nPatterns,1);
meanRelAmps = zeros(nPatterns, size(details.actualRelAmps{1},1));
for ii = 1:nPatterns
    if  sum(any(details.sAmps{ii},2) ~= 0) == 1 && ~details.sAloneBin(ii) %electrode pairs
        pairFlags(ii) = true;
    end
    meanRelAmps(ii,:) = mean(details.actualRelAmps{ii},2);
end

tripletCombs = struct([]);
ampDiffs = [];
normAmpDiffs = [];
sAloneInd = find(details.sAloneBin)';



for ii = 1:nPatterns
    if sum(any(details.sAmps{ii},2) ~= 0) == 2 % electrode triplets
        sElecsActive = find(any(details.sAmps{ii},2));
        tripletCombs(end+1).tripInd = ii; %#ok<AGROW>
        for jj = 1:2
            relAmp = meanRelAmps(ii,sElecsActive(jj));
            %find paired electrode pattern with mean relative amplitude
            %within 5% of mean relative amplitude used in triplet
            pairPattern =  find((abs(meanRelAmps(:,sElecsActive(jj))-relAmp) < abs(0.05*relAmp)).*pairFlags);
            tripletCombs(end).pairInds(jj) = pairPattern;
            ampDiffs = [ampDiffs abs(meanRelAmps(pairPattern, sElecsActive(jj))-relAmp)];
            normAmpDiffs = [normAmpDiffs abs(meanRelAmps(pairPattern, sElecsActive(jj))-relAmp)/abs(relAmp)];
            
            %find corresponding secondary alone pattern and retrieve
            %threshold
            found = false;
            for kk = sAloneInd
                if any(details.sAmps{kk}(sElecsActive(jj),:)) &&...
                    sign(details.sAmps{ii}(sElecsActive(jj),1)) == sign(details.sAmps{kk}(sElecsActive(jj),1))
                    tripletCombs(end).sAloneThreshes(jj) = thresholds(kk);
                    if isfield(details, 'threshMins')
                        tripletCombs(end).sAloneThreshMins(jj) = details.threshMins(kk);
                    end
                    if found
                        error('found 2 matching s-alone stimuli--this shouldn''t ever happen!')
                    end
                    found = true;
                end
            end
            if ~found
                error('didn''t find any matching s-alone stimuli--this shouldn''t ever happen!')
            end
        end
        
        %flag any triplets containing a component pair that wasn't used in
        %linear fits of electrode pair pattern thresholds (because of infinite threshold,
        %specific exclusion in cell_list_cluster_stim, or rounding of some secondary amplitudes to 0) 
        if details.unfitPairs(tripletCombs(end).pairInds(1)) || details.unfitPairs(tripletCombs(end).pairInds(2));
            tripletCombs(end).unfitPair = true;
        else
            tripletCombs(end).unfitPair = false;
        end
        
        %flag any triplets containing a component pair that doesn't fall
        %into positive primary region
        tripletCombs(end).inPrimReg = true;
        for jj = 1:2
            mBound = details.mValidBound(:,sElecsActive(jj));
            relAmp = meanRelAmps(ii,sElecsActive(jj));
            if (mBound(1) && relAmp > 0 && (1/relAmp) < mBound(1)) || (mBound(2) && relAmp < 0 && abs(1/relAmp) < abs(mBound(2)))
                tripletCombs(end).inPrimReg = false;
                %disp([num2str(meanRelAmps(ii,:)) ' (' num2str(relAmp) ')'])
            end
        end
    end
end


if ~isempty(estPrimThresh)
    primThresh = estPrimThresh;
elseif useFitPAloneThresh %use mean of the y-intercepts of the
    primThresh = details.tp;
    %primThresh = mean(details.yIntercepts);
else %use actual measured primary-alone threshold
    pAloneInd = find(details.pAloneBin);
    pAloneAmps = [details.pAmps{pAloneInd(1)}(1) details.pAmps{pAloneInd(2)}(1)];
    primThresh = thresholds(pAloneInd(pAloneAmps<0));
end

disp(['relative amplitude difference between pairs and triplets: '...
    num2str(min(ampDiffs)) '-' num2str(max(ampDiffs)) ' (mean = ' num2str(mean(ampDiffs)) ')'])

disp(['fractional relative amplitude difference between pairs and triplets: '...
    num2str(100*min(normAmpDiffs)) '%-' num2str(100*max(normAmpDiffs)) '% (mean = ' num2str(100*mean(normAmpDiffs)) '%)'])

%%

nCombs = length(tripletCombs);
xToPlot = zeros(nCombs, 1);
yToPlot = zeros(nCombs, 1);
zToPlot = zeros(nCombs, 1);
predThresh = zeros(nCombs,1);
shiftExpected = zeros(nCombs, 1);

for ii = 1:nCombs
    
    if ~fullModel %estimated triplet thresholds shift is based on measured shifts of component pairs
        xToPlot(ii) = thresholds(tripletCombs(ii).pairInds(1)) - primThresh;
        yToPlot(ii) = thresholds(tripletCombs(ii).pairInds(2)) - primThresh;
        zToPlot(ii) = thresholds(tripletCombs(ii).tripInd) - primThresh;
    else %estimated triplet thresholds shift is based on slopes of lines fit to component pairs
        %mean relative amplitudes on component pairs
        sInd1 = find(meanRelAmps(tripletCombs(ii).pairInds(1),:));
        sInd2 = find(meanRelAmps(tripletCombs(ii).pairInds(2),:));
%         relAmp1 = meanRelAmps(tripletCombs(ii).pairInds(1),sInd1);
%         relAmp2 = meanRelAmps(tripletCombs(ii).pairInds(2),sInd2);
        relAmp1 = meanRelAmps(tripletCombs(ii).tripInd,sInd1);
        relAmp2 = meanRelAmps(tripletCombs(ii).tripInd,sInd2);
        
        
        %find I_i at intercept between fit line of model (yIntercept-lambda_i*I_i) and line representing current combinations
        %used in the pattern (1/relAmps)*I_i
        yInt1 = details.tp; lambda1 = details.lambdas(sInd1);
        I_1 = yInt1/((1/relAmp1)+lambda1);
        
        yInt2 = details.tp; lambda2 = details.lambdas(sInd2);
        I_2 = yInt2/((1/relAmp2)+lambda2);      
        
        %shift = -1*lambda_i*I_i
        xToPlot(ii) = -lambda1*I_1;
        yToPlot(ii) = -lambda2*I_2;

        %temporary - used to confirm accuracy of other code
        %predThresh(ii) = xToPlot(ii)+yToPlot(ii)+primThresh;
        
        %if lines intercept at a point corresponding to I_0 < 0, set to
        %infinity (this should only happen for large abs(relAmps))
        if I_1/relAmp1 < 0
            xToPlot(ii) = inf;
            disp('model line doesn''t intersect line representing current ratio')
            keyboard
        end
        if I_2/relAmp2 < 0
            yToPlot(ii) = inf;
            disp('model line doesn''t intersect line representing current ratio')
            keyboard
        end
        zToPlot(ii) = thresholds(tripletCombs(ii).tripInd) - primThresh;
        
        
                
        %%% corrected version %%%
        r1 = relAmp1;
        r2 = relAmp2;
        
        tp = details.tp;
        
        shiftExpected(ii) = (tp/(1+lambda1*r1+lambda2*r2)) - tp;
        
        %%%
        
    end
    
    
    
    %make x-axis always denote secondary electrode that has larger (positive) effect (since they are otherwise
    %arbitrarily assigned)
    if xToPlot(ii) < yToPlot(ii)
        tempX = yToPlot(ii);
        yToPlot(ii) = xToPlot(ii);
        xToPlot(ii) = tempX;
    end
end

%convert thresholds shifts from absolute current amplitude to percent of p-alone amplitude
xToPlot = 100*xToPlot/primThresh;
yToPlot = 100*yToPlot/primThresh;
zToPlot = 100*zToPlot/primThresh;

shiftExpected = 100*shiftExpected/primThresh;

%separate out any points in which one of the coordinates is infinity
xToPlotInf = []; yToPlotInf = []; zToPlotInf = []; tripletCombsInf = []; %for valid x and y values but infinite z value
xToPlotInvalid = []; yToPlotInvalid = []; zToPlotInvalid = []; tripletCombsInvalid = []; %if x or y value is infinite

for ii = nCombs:-1:1
    if any(isinf([xToPlot(ii) yToPlot(ii) zToPlot(ii)]))
        disp(['one of the following pattern indeces had infinite threshold: ',...
            10, num2str(tripletCombs(ii).pairInds(1)) ', ' num2str(tripletCombs(ii).pairInds(2)) ', '...
            num2str(tripletCombs(ii).tripInd)])
        if isinf(xToPlot(ii)) || isinf(yToPlot(ii))
            disp(['*****one of the following component pair pattern indeces has infinite threshold: ',...
                10, num2str(tripletCombs(ii).pairInds(1)) ', ' num2str(tripletCombs(ii).pairInds(2)) ': excluding from plot!*****'])
            xToPlotInvalid = [xToPlot(ii) xToPlotInvalid]; %#ok<AGROW>
            yToPlotInvalid = [yToPlot(ii) yToPlotInvalid]; %#ok<AGROW>
            zToPlotInvalid = [zToPlot(ii) zToPlotInvalid]; %#ok<AGROW>
            tripletCombsInvalid = [tripletCombs(ii) tripletCombsInvalid]; %#ok<AGROW>
        else %only zToPlot is infinite
            disp(['measured triplet threshold (pattern index ' num2str(tripletCombs(ii).tripInd) ') is set to infinity'])
            
            xToPlotInf = [xToPlot(ii) xToPlotInf]; %#ok<AGROW>
            yToPlotInf = [yToPlot(ii) yToPlotInf]; %#ok<AGROW>
            zToPlotInf = [zToPlot(ii) zToPlotInf]; %#ok<AGROW>
            tripletCombsInf = [tripletCombs(ii) tripletCombsInf]; %#ok<AGROW>
        end
        xToPlot(ii) = [];
        yToPlot(ii) = [];
        zToPlot(ii) = [];
        tripletCombs(ii) = [];
        
        shiftExpected(ii) = [];
    end
end

nCombsUsed = length(xToPlot);


%% 2-D figure with colormapping: may not work (several code changes since
% last use)
if 0
    xRange = max(xToPlot) - min(xToPlot);
    yRange = max(yToPlot) - min(yToPlot);
    xMed = mean([max(xToPlot) min(xToPlot)]);
    yMed = mean([max(yToPlot) min(yToPlot)]);
    
    halfPlotLength = max(xRange, yRange)*0.55;

    xPlotLimits = [-halfPlotLength halfPlotLength] + xMed;
    yPlotLimits = [-halfPlotLength halfPlotLength] + yMed;
    
    %least-squares fitting of plane to data,
    A = ones(nCombsUsed, 3);
    A(:,1) = xToPlot;
    A(:,2) = yToPlot;
    b(:,1) = zToPlot;
    
    %best-fit plane: c(1)*x + c(2)*y + c(3) = z
    c = A\b;
    
    xGrid = linspace(xPlotLimits(1), xPlotLimits(2), 10);
    yGrid = linspace(yPlotLimits(1), yPlotLimits(2), 10);
    zVals = zeros(length(yGrid), length(xGrid));
    for i = 1:length(xGrid)
        for j = 1:length(yGrid)
            zVals(j,i) = c(1)*xGrid(i) + c(2)*yGrid(j) + c(3);
        end
    end
    
    zMax = max(zToPlot);
    zMin = min(zToPlot);
    
    colors = jet(501);
    figure('position', [100 100 500 500])
    hold on
    
    %zContourLevels = linspace(zMin + (zMax-zMin)/20, zMax - (zMax-zMin)/20, 9);
    zContourLevels = linspace(zMin, zMax, 6);
    
    for ii = 1:length(zContourLevels)
        colorToPlot = floor(500*(zContourLevels(ii)-zMin)/(zMax-zMin))+1;
        contour(xGrid, yGrid, zVals, [zContourLevels(ii) zContourLevels(ii)], 'color', colors(colorToPlot,:))
    end
    
    for ii = 1:nCombsUsed
        colorToPlot = floor(500*(zToPlot(ii)-zMin)/(zMax-zMin))+1;
        plot(xToPlot(ii), yToPlot(ii), 'o', 'markerFaceColor', colors(colorToPlot, :), 'markerEdgeColor', 'k')
    end
    
    %plot points corresponding to infinite z values (threshold above measured range)
    for ii = 1:length(xToPlotInf)
        if zToPlotInf(ii) > 0
            plot(xToPlotInf(ii), yToPlotInf(ii), 'o', 'markerFaceColor', 'w', 'markerEdgeColor', 'k')
            plot(xToPlotInf(ii), yToPlotInf(ii), 'k+')
        else %z is negative infinity
            plot(xToPlotInf(ii), yToPlotInf(ii), 'o', 'markerFaceColor', 'w', 'markerEdgeColor', 'k')
            text(xToPlotInf(ii), yToPlotInf(ii), '-')
        end
    end
    
    
    %contour(xGrid, yGrid, zVals, 10)
    hold off
    xlabel('percent shift in threshold due to secondary electrode 1')
    ylabel('percent shift in threshold due to secondary electrode 2')
    title('shifts in threshold from estimated primary-alone threshold')
    %set(gca, 'xLim', [min(xToPlot) max(xToPlot)], 'yLim', [min(yToPlot) max(yToPlot)])
    
    tickDiv = 10;
    xTickLims = [tickDiv*ceil(xPlotLimits(1)/tickDiv) tickDiv*floor(xPlotLimits(2)/tickDiv)];
    yTickLims = [tickDiv*ceil(yPlotLimits(1)/tickDiv) tickDiv*floor(yPlotLimits(2)/tickDiv)];
    xTickLocs = xTickLims(1):tickDiv:xTickLims(2);
    yTickLocs = yTickLims(1):tickDiv:yTickLims(2);
    set(gca, 'xLim', [xPlotLimits(1) xPlotLimits(2)], 'yLim', [yPlotLimits(1) yPlotLimits(2)], 'box', 'on',...
        'xTick', xTickLocs, 'yTick', yTickLocs)
    axis equal
    
    
    %plot vertical color scale bar
    colors = jet(100);
    figure('position', [100 100 100 500])
    hold on
    nColors = size(colors, 1);
    xPoly = [0 1 1 0];
    for i = 1:nColors
        yPoly = [(i-1)/nColors (i-1)/nColors (i+0.05)/nColors (i+0.05)/nColors];
        fill(xPoly, yPoly, colors(i,:), 'EdgeColor', colors(i,:), 'LineWidth', 0.001)
    end
    
    if zMax-zMin < 100
        colorBarLabels = -100:10:100;
    elseif zMax-zMin < 250
        colorBarLabels = -300:20:300;
    else
        colorBarLabels = -500:50:500;
    end
    for i = 1:length(colorBarLabels)
        if colorBarLabels(i) > zMin && colorBarLabels(i) < zMax
            y = (colorBarLabels(i)-zMin)/(zMax-zMin);
            plot([1 1.1], [y y], 'k-')
            text(1.2, y, [num2str(colorBarLabels(i)), '%'])
        end
    end
    set(gca, 'ylim', [0 1], 'xlim', [0 2])
    axis off
    hold off
end

%% expected vs. observed plot

if isempty(obsVExpAxes)
    figure('position', [100 100 300 300])
    axes('position', [0.2 0.2 0.7 0.7])
else
    axes(obsVExpAxes)
end
hold on

if ~isempty(obsVExpColor)
    plotColor = obsVExpColor;
else
    plotColor = [0 0 0];
end


unfitCompPairs = [tripletCombs.unfitPair];
compPairOutsidePrimRegion = ~[tripletCombs.inPrimReg];


%flags triplets patterns that are excluded for a
%reason other than just being outside positive primary region
badCompPairBin = false(nPatterns,1);
badCompPairBin([tripletCombs(unfitCompPairs).tripInd]) = true; %indexed by Pattern (vs. tripletCombs indexing)

%exclude all patterns corresponding to component pairs that are
%specifically flagged as bad in cell_list_cluster_stim, pairs
%that fall outside positive primary region, or pairs that weren't used in
%pair model fits (e.g. because some secondary amplitudes rounded to 0 or because infinite thresh)

excludeFromFits = unfitCompPairs | compPairOutsidePrimRegion;

unfitTripBin = false(nPatterns,1);
unfitTripBin([tripletCombs(excludeFromFits).tripInd]) = true; %indexed by pattern (vs. tripletCombs indexing)

%plot symbols depict component pair validity
toPlotO = ~excludeFromFits; %component pairs within valid zone and none of 3 patterns flagged as bad
toPlotX = unfitCompPairs | compPairOutsidePrimRegion; %at least one component pair wasn't used in fitting in positive primary fitting region
toPlotStar = excludeFromFits & ~(unfitCompPairs | compPairOutsidePrimRegion); %both component were used in fitting in positive primary fitting region but 
% one of the 3 patterns specifically flagged as bad in cell_list_cluster_stim.m

if ~plotUnfitted
    toPlotX = false(size(toPlotX)); %should all be included in excludeFromFits    
    toPlotStar = false(size(toPlotStar)); %should all be included in excludeFromFits    
end


toPlotSum = zeros(size(toPlotO)); %to check for overlap/gaps between categories

toPlot = toPlotO; toPlotSum = toPlotSum + toPlot;
plot(shiftExpected(toPlot), zToPlot(toPlot), '.', 'MarkerEdgeColor', plotColor)
%plot(xToPlot(toPlot)+yToPlot(toPlot), zToPlot(toPlot), '.', 'MarkerEdgeColor', plotColor)

toPlot = toPlotX; toPlotSum = toPlotSum + toPlot;
plot(shiftExpected(toPlot), zToPlot(toPlot), 'x', 'MarkerEdgeColor', plotColor)
%plot(xToPlot(toPlot)+yToPlot(toPlot), zToPlot(toPlot), 'x', 'MarkerEdgeColor', plotColor)

toPlot = toPlotStar; toPlotSum = toPlotSum + toPlot;
plot(shiftExpected(toPlot), zToPlot(toPlot), '*', 'MarkerEdgeColor', plotColor)
%plot(xToPlot(toPlot)+yToPlot(toPlot), zToPlot(toPlot), '*', 'MarkerEdgeColor', plotColor)


%check for plot category overlap/gaps - if so there is an error in the code!
if ~plotUnfitted
    if ~all((toPlotSum + excludeFromFits) == 1);
        error('gaps or overlaps between plot categories')
    end        
elseif ~all(toPlotSum == 1);
    error('gaps or overlaps between plot categories')
end

%least-squares fitting of line to data,
A = ones(nCombsUsed, 2);
A(excludeFromFits,:) = [];

%A(:,1) = xToPlot(~excludeFromFits)+yToPlot(~excludeFromFits);
A(:,1) = shiftExpected(~excludeFromFits);
b(:,1) = zToPlot(~excludeFromFits);

c = A\b;


%xDataLims = [min(xToPlot(toPlotSum~=0)+yToPlot(toPlotSum~=0)) max(xToPlot(toPlotSum~=0)+yToPlot(toPlotSum~=0))];
xDataLims = [min(shiftExpected(toPlotSum~=0)) max(shiftExpected(toPlotSum~=0))];
yDataLims = [min(zToPlot(toPlotSum~=0)) max(zToPlot(toPlotSum~=0))];

xDataMed = mean(xDataLims);
yDataMed = mean(yDataLims);

%plot the best-fit line and the model prediction line
plot(xDataLims, xDataLims*c(1)+c(2), '-', 'color', plotColor)
plot([min(zToPlot) max(zToPlot)], [min(zToPlot) max(zToPlot)], 'k--')

% calculate R-squared (where model is observed = predicted from additivity, vs. observed = predicted by best linear fit)
%pred = xToPlot(~excludeFromFits)+yToPlot(~excludeFromFits);
pred = shiftExpected(~excludeFromFits);
obs = zToPlot(~excludeFromFits);
r_sq = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);

%plots points corresponding to infinite z values (threshold above measured
%range)
if plotUnfitted
    warned = false;
    for ii = 1:length(xToPlotInf)
        if ~warned
            warning('markers of infinite z values do not reflect categories used for noninfinite values')
            warned = true;
        end
        if zToPlotInf(ii)>0
            plot(xToPlotInf(ii) + yToPlotInf(ii), max(zToPlot), 'k^')
        else %z = -inf
            plot(xToPlotInf(ii) + yToPlotInf(ii), min(zToPlot), 'kv')
        end
    end
end

hold off
xlabel('expected percent threshold shift')
ylabel('observed percent threshold shift')

tickDiv = 20;

halfPlotLength = max([xDataLims(2)-xDataLims(1) yDataLims(2)-yDataLims(1)])*0.55;
xPlotLimits = [-halfPlotLength halfPlotLength] + xDataMed;
yPlotLimits = [-halfPlotLength halfPlotLength] + yDataMed;


xTickLims = [tickDiv*ceil(xPlotLimits(1)/tickDiv) tickDiv*floor(xPlotLimits(2)/tickDiv)];
yTickLims = [tickDiv*ceil(yPlotLimits(1)/tickDiv) tickDiv*floor(yPlotLimits(2)/tickDiv)];
xTickLocs = xTickLims(1):tickDiv:xTickLims(2);
yTickLocs = yTickLims(1):tickDiv:yTickLims(2);

axis equal
set(gca, 'xLim', xPlotLimits, 'yLim', yPlotLimits, 'xTick', xTickLocs, 'yTick', yTickLocs, 'box', 'on')


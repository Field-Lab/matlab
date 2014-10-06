function [thresholds threshStds estimatedPrimThresh details] =...
    generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID, relAmps, varargin)
%
%
% arguments:
%   pathToData  should be path to elecResp files, which must be named in the format: 'elecResp_n
%               [neuronID]  _p [pattern number]'
%   patternNos  vector of integers corresponding to pattern numbers in elecResp names
%   pElec       single integer specifying primary electrode
%   neuronID    
%
%   relAmps     amplitudes on secondary electrodes, relative to primary 
%               ***should be listed from smallest value (e.g. -1) to highest value (e.g. 1) and NOT
%               include the value 0
%   varargin    see parameter list in m-file
%
% outputs:
%   thresholds  vector of threshold values with indices corresponding to patternNos
%   threshStds  vector of threshold standard deviations based on bootstrapping, with indices
%               corresponding to patternNos
%   estimatedPrimThresh
%               estimated threshold to primary electrode alone (cathodal) stimulation, based on
%               linear regression of primary+single secondary thresholds
%   details
%               .sAmps     cell array of vectors of lowest amplitude on secondary electrode for
%                          given pattern: details.sAmps{pattern number}(index of secondary
%                          electrode, corresponding to sElecs)
%               .actualRelAmps
%                          organized as .sAmps, but divided by corresponding .pAmps (secondary
%                          electrode amplitude relative to primary electrode amplitude)
%                          
%               .pAmps     vector of lowest amplitudes on primary electrode, with indices
%                          corresponding to patternNos
%               .sElecs    vector of secondary electrodes
%               .pElec     primary electrode
%               .pAloneBin
%                          binary vector with 'true' for pattern number corresponding to
%                          primary-electrode-alone stimulation
%               .curveWidths
%                          vector of values indicating second parameter defining cumulative gaussian
%                          fit, where the standard deviation of the guassian = threshold/width
%               .threshShiftSlopes
%                          will be phased out
%               .secondaryCombsIncluded
%                          cell array of 3-element vectors specifying, for each non-excluded pattern
%                          with a triplet of electrodes, (1) pattern number corresponding to
%                          triplet, (2,3) pattern numbers corresponding to component pairs of
%                          triplet
%                                e.g.
%                                     pattern 50 corresponds with stim triplet with secondary elec 1 relative
%                                        amp 2x, and secondary electrode 2 relative amp 0.5x
%                                     pattern 3 corresponds with stim pair with secondary elec 1
%                                        relative amp 2x
%                                     pattern 7 corresponds with stim pair with secondary elec 2
%                                        relative amp 0.5x
%                                     .secondaryCombsIncluded{x} = [50 3 7]
%
% CHANGES
% 2010-04-07
% changed calculation of estimated primary-alone threshold to reflect better linear fit in terms of
% absolute (vs relative) amplitude on secondary
%
% 2011-11-30
% due to variations in relative secondary amplitude with movie number, mean
% is calculated over all movies to determine values in details.actualRelAmps and
% standard deviations are stored in details.actualRelAmpsSD
% * pAmps and sAmps still reflect amplitudes found in first movie




p = inputParser;

p.addRequired('pathToData', @ischar)
p.addRequired('patternNos', @isnumeric)
p.addRequired('pElec', @isnumeric)
p.addRequired('neuronID', @isnumeric)
p.addRequired('relAmps', @isnumeric)

p.addParamValue('bootstrapReps', 100, @isnumeric)
p.addParamValue('plotLogScale', 0, @(x)any(x==[0 1])) %% plots thresholds and curves in log scale
p.addParamValue('constrainSlopes', 0, @(x)any(x==[0 1]))
p.addParamValue('redoConstrainedFitting', 0, @(x)any(x==[0 1]))
p.addParamValue('threshPlotLims', [0 1.5], @isnumeric)
p.addParamValue('curvePlotXLims', [0 1.5], @isnumeric)
p.addParamValue('ampMatchTol', 0.005, @isnumeric)
p.addParamValue('plotEstimatedPrimThresh', 0, @(x)any(x==[0 1]))
p.addParamValue('excludePatterns', [], @isnumeric) %excludes patterns from constrained erf slope analysis, erf curves not plotted in summary plot or isoresponse contour plots, and triplet patterns containing secondary electrodes amplitudes are marked for exclusion from later analysis
p.addParamValue('recalcAll', 0, @(x)any(x==[0 1]))
p.addParamValue('insetYLim', [0 1], @isnumeric)
p.addParamValue('erfStartParams', [], @isnumeric)
p.addParamValue('plotInconsistentFits', false, @islogical)

p.parse(pathToData, patternNos, pElec, neuronID, relAmps, varargin{:})

nBootstrapReps = p.Results.bootstrapReps;
plotLogScale = p.Results.plotLogScale;
%constrainSlopes = p.Results.constrainSlopes;
redoConstrainedFitting = p.Results.redoConstrainedFitting;
%threshPlotLims = p.Results.threshPlotLims;
curvePlotXLims = p.Results.curvePlotXLims;
ampMatchTol = p.Results.ampMatchTol;
%plotEstimatedPrimThresh = p.Results.plotEstimatedPrimThresh;
excludePatterns = p.Results.excludePatterns;
recalcAll = p.Results.recalcAll;
insetYLim = p.Results.insetYLim;
erfStartParams = p.Results.erfStartParams;
plotInconsistentFits = p.Results.plotInconsistentFits; %gets passed to checkForUnfinishedAnalysis

cd(pathToData)

nPatterns = length(patternNos);
nRelAmps = length(relAmps);

sElecs = getCluster(pElec);
sElecs = sElecs(2:end); %first is primary electrode

elecCombs =              {};
combIDs =                zeros(nPatterns, 1);
elecs =                  cell(nPatterns, 1);
baseAmps =               cell(nPatterns, 1);
erfParams =              cell(nPatterns, 1);
erfErrors =              zeros(nPatterns, 1);
thresholds =             zeros(nPatterns, 1);
threshStds =             zeros(nPatterns, 1);
curveWidths =            zeros(nPatterns, 1);
pAmps =                  zeros(nPatterns, 1);
sAmps =                  cell(nPatterns, 1); % each cell is a vector corresponding to sElecs
currentAmpLims =         cell(nPatterns, 1);
allData =                cell(nPatterns, 1);
actualRelAmps =          cell(nPatterns, 1);
actualRelAmpsSD =        cell(nPatterns, 1);
pAloneBin =              zeros(nPatterns, 1);
sAloneBin =              zeros(nPatterns, 1);


for i = 1:nPatterns
    
    disp(num2str(i))
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
    elecResp = temp.elecResp;
    
    elecResp.stimInfo.pElec = pElec;
    
    nMovies = length(elecResp.stimInfo.movieNos);
    neuronID = elecResp.cells.main;
    pathToEi = elecResp.names.rrs_ei_path;
    

    analyzedMovies = ones(nMovies, 1);
    for j = 1:nMovies
        if isempty(elecResp.analysis.type{j})
            analyzedMovies(j) = 0;
        end
    end
    
    
    if ~isempty(erfStartParams)
        elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll,...
            'erfStartParams', erfStartParams, 'plotInconsistentFits', plotInconsistentFits);
    else        
        elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll,...
           'plotInconsistentFits', plotInconsistentFits);
    end
    
    save(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))], 'elecResp')
    
    
    %extracts info necessary to plot

    %allData only used for slider plot of secondary combination curves
    allData{i} = zeros(2, nMovies);
    allData{i}(1,:) = abs(elecResp.stimInfo.stimAmps);
    allData{i}(2,:) = elecResp.analysis.successRates;
    for j = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{j})
            allData{i}(:,j) = [];
        end
    end
    
%     if constrainSlopes %currently broken (only works for log-scale erf)
%         erfParams{i} = elecResp.analysis.constrainedSlopeLogErfParams;
%         thresholds(i) = elecResp.analysis.constrainedSlopeLogThresh;
%         erfErrors(i) = elecResp.analysis.constrainedSlopeLogErfErr;
%         threshStds(i) = elecResp.analysis.constrainedSlopeLogThreshStd;
%     else
        erfParams{i} = elecResp.analysis.erfParams; %individually  fit erfs
        thresholds(i) = elecResp.analysis.threshold; %based on individually fit erfs
        threshStds(i) = elecResp.analysis.threshStd;
        erfErrors(i) = elecResp.analysis.erfErr;
        curveWidths(i) = -sqrt(2)*elecResp.analysis.erfParams(2);
%     end

    % determines plot x-limits for each individual response curve
    analyzedMoviesI = find(analyzedMovies);    
    currentAmpLims{i} = [elecResp.stimInfo.stimAmps(analyzedMoviesI(1))...
        elecResp.stimInfo.stimAmps(analyzedMoviesI(end))];
    
    
    % sorts out stimulus information
    elecs{i} = elecResp.stimInfo.electrodes;
    baseAmps{i} = getStimAmps(elecResp.names.data_path, elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(1));
    
    %checks to makes sure that ratio is consistent
    %relAmpsTemp = baseAmps{i}/baseAmps{i}(1);
    
    % amplitude on primary
    if any(elecs{i} == pElec)
        pAmps(i) = baseAmps{i}(elecs{i} == pElec);
    else
        pAmps(i) = 0;
    end
        
    %amplitude on secondary(ies)
    sAmps{i} = zeros(length(sElecs), 1);
    actualRelAmps{i} = zeros(length(sElecs), 1); %sometimes deviate slightly from ideal (e.g. 0.5, 1, 2)

    alreadyWarned = false;
    for k = 1:length(elecs{i})
        if any(elecs{i} == pElec)
            
            sInd = find(sElecs==elecs{i}(k));
            if ~isempty(sInd) %if current electrode is a secondary electrode
                sAmps{i}(sInd) = baseAmps{i}(k);
                %actualRelAmps{i}(j) = sAmps{i}(j)/pAmps(i);
                
                relAmpsAllMov = zeros(1,nMovies);
                for jj = 1:nMovies
                    [tempAmps tempElecs] = getStimAmps(elecResp.names.data_path, elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(jj));
                    if any(tempElecs == sElecs(sInd)) %secondary amplitude sometimes rounds to zero
                        relAmpsAllMov(jj) = tempAmps(tempElecs == sElecs(sInd))/tempAmps(tempElecs == pElec);
                        
%                         if ~all(tempAmps == relAmpsTemp)
%                             warning(['relative amplitudes are not consistent throughout all movies of pattern ' num2str(patternNos(i))])
%                             break
%                         end
                    elseif ~alreadyWarned %secondary amplitude rounded to zero
                        disp(['warning: at least one secondary amplitude in pattern ' num2str(patternNos(i)) ' rounded to zero'])
                        alreadyWarned = true;
                    end
                end
                
                actualRelAmps{i}(sInd) = mean(relAmpsAllMov);
                actualRelAmpsSD{i}(sInd) = std(relAmpsAllMov);
            end
        else
            actualRelAmps{i}(sElecs==elecs{i}(k)) = inf;
            sAmps{i}(sElecs==elecs{i}(k)) = baseAmps{i}(k);
        end
    end
    
    if pAmps(i) && ~any(sAmps{i})
        pAloneBin(i) = 1;
    end
    
    if any(sAmps{i}) && ~pAmps(i)
        sAloneBin(i) = 1;
    end
    
    % adds combination of electrodes to elecCombs list if not already in list
    foundMatch = 0;
    for j = 1:length(elecCombs)
        if length(elecs{i}) == length(elecCombs{j})
            if all(elecs{i} == elecCombs{j}) %%should fix so that orders don't have to match
                foundMatch = 1;
                combIDs(i) = j; %maps pattern index to elecComb index
                break
            end
        end
    end
    
    if ~foundMatch
        elecCombs{length(elecCombs)+1} = elecs{i}; %#ok<AGROW>
        combIDs(i) = length(elecCombs); %index = last in current elecCombs
    end
    
end

details.sAmps = sAmps;
details.actualRelAmps = actualRelAmps;
details.pAmps = pAmps;
details.sElecs = sElecs;
details.pElec = pElec;
details.pAloneBin = pAloneBin;
details.sAloneBin = sAloneBin;
details.curveWidths = curveWidths;

%% fitting with constrained slopes: needs to be changed for compliance with linear erf fitting

if redoConstrainedFitting
    elecRespArray = fitErfsConstrainedSlopeWrapper(patternNos, neuronID, pathToData, nBootstrapReps, excludePatterns);
    
    pIndex = 0;
    for i = 1:nPatterns
        if ~any(excludePatterns == patternNos(i))
            pIndex = pIndex + 1;
            elecResp = elecRespArray{pIndex};
            save(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))], 'elecResp')
        end
    end
end

%% figure

figure('position', [100 100 1000 800], 'color', 'white')

centerAxes = axes('position', [0.3 0.3 0.4 0.4]);


%starts at top and goes clockwise
% outerAxes{1} = axes('position', [0.375 0.75 0.25 0.2]);
% outerAxes{2} = axes('position', [0.7 0.55 0.25 0.2]);
% outerAxes{3} = axes('position', [0.7 0.25 0.25 0.2]);
% outerAxes{4} = axes('position', [0.375 0.05 0.25 0.2]);
% outerAxes{5} = axes('position', [0.1 0.25 0.25 0.2]);
% outerAxes{6} = axes('position', [0.1 0.55 0.25 0.2]);

outerAxes{1} = axes('position', [0.7 0.7 0.25 0.2]);
outerAxes{2} = axes('position', [0.7 0.4 0.25 0.2]);
outerAxes{3} = axes('position', [0.7 0.1 0.25 0.2]);

outerAxes{4} = axes('position', [0.1 0.1 0.25 0.2]);
outerAxes{5} = axes('position', [0.1 0.4 0.25 0.2]);
outerAxes{6} = axes('position', [0.1 0.7 0.25 0.2]);

pInd = strfind(elecResp.names.rrs_short_name, 'p');
text(0, 2, elecResp.names.rrs_short_name(1:pInd(end)-2), 'units', 'normalized', 'interpreter', 'none')
meanSqError = mean(erfErrors.^2);
text(0, 1.9, ['mean squared error: ' num2str(meanSqError)], 'units', 'normalized', 'interpreter', 'none')

if nRelAmps == 2
    colorScheme(1, :) = [1 0 0]; % red
    colorScheme(2, :) = [0 0 1];
elseif nRelAmps == 6
    colorScheme(1, :) = [1 0.5 0.5]; %pink
    colorScheme(2, :) = [1 0 0]; % red
    colorScheme(3, :) = [0.75 0 0]; %dark red
    colorScheme(4, :) = [0 0 0.75]; %dark blue
    colorScheme(5, :) = [0 0 1]; %blue
    colorScheme(6, :) = [0.5 0.5 1]; %light blue
elseif nRelAmps == 8
    colorScheme(1, :) = [1 0.5 0.5]; %pink
    colorScheme(2, :) = [1 0 0]; % red
    colorScheme(3, :) = [0.75 0 0]; %dark red
    colorScheme(4, :) = [0.5 0 0]; %darker red
    colorScheme(5, :) = [0 0 0.5]; %darker blue
    colorScheme(6, :) = [0 0 0.75]; %dark blue
    colorScheme(7, :) = [0 0 1]; %blue
    colorScheme(8, :) = [0.5 0.5 1]; %light blue
else
    disp('Only 6 or 8 relative amplitudes are allowed!!!')
    return
end


% legend

axes('position', [0.43 0.72 0.15 0.18])
hold on

text(1, 9.8, 'secondary amplitude,')
text(1, 9, 'relative to primary')

markerPositions = linspace(8, 2, nRelAmps+1);
plot(1, markerPositions(1), 's', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0])
text(2, markerPositions(1), num2str(0))

for i = 1:nRelAmps
    plot(1, markerPositions(i+1), 's', 'MarkerFaceColor', colorScheme(i,:), 'MarkerEdgeColor', colorScheme(i,:))
    %text(2, markerPositions(i+1), ['~' num2str(relAmps(i))])
    text(2, markerPositions(i+1), num2str(relAmps(i)))
end

set(gca, 'xlim', [0 10], 'ylim', [0 10])
axis off
hold off

%% plots ei and determines which secondary electrode corresponds with which axis position

%plotEi61(pathToEi, neuronID, 'axesH', centerAxes, 'markElecs', pElec)
[xCoords yCoords] = getElectrodeCoords61();


% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
clear eiFile

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for j = 1:64
    if ~(j==9||j==25||j==57)
        eiAmps(j) = max(max(abs(ei(1,j+1,:))));
    end
end


eiAmps = eiAmps/max(eiAmps);
hex_contour(xCoords, yCoords, eiAmps, 8, 'fig_or_axes', centerAxes, 'contourSpacing', 'linear', 'plotCoords', false);
hold on

%primary electrode
plot(xCoords(pElec), yCoords(pElec), 'k*')

%other electrodes
for i = 1:64
    if ~any([9 25 57 pElec] == i)
        plot(xCoords(i), yCoords(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
    end
end


%array outline
plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')

set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
axis equal
axis off
hold off




[xCoords yCoords] = getElectrodeCoords61();

pPos = [xCoords(pElec) yCoords(pElec)];
secondaryAxes = zeros(6, 1); %maps sElecs onto axes position ids, starting at top and going clockwise
for i = 1:6
    sPos = [xCoords(sElecs(i)) yCoords(sElecs(i))];
    if sPos(1) == pPos(1) && sPos(2) > pPos(2)
        secondaryAxes(i) = 1;
    elseif sPos(1) > pPos(1) && sPos(2) > pPos(2)
        secondaryAxes(i) = 2;
    elseif sPos(1) > pPos(1) && sPos(2) < pPos(2)
        secondaryAxes(i) = 3;
    elseif sPos(1) == pPos(1) && sPos(2) < pPos(2)
        secondaryAxes(i) = 4;
    elseif sPos(1) < pPos(1) && sPos(2) < pPos(2)
        secondaryAxes(i) = 5;
    else
        secondaryAxes(i) = 6;
    end
end


threshCurves = cell(6,1);
threshCurveStds = cell(6,1);
threshCurveActualRelAmps = cell(6,1);
for i = 1:6
    threshCurves{i} = inf*ones(nRelAmps, 1);
    threshCurveStds{i} = inf*ones(nRelAmps, 1);
    threshCurveActualRelAmps{i} = inf*ones(nRelAmps, 1);
end


for i = 1:nPatterns
    if ~any(excludePatterns == patternNos(i))
        
        xProj = abs(currentAmpLims{i}(1)):0.01:abs(currentAmpLims{i}(2));
        projection = 0.5 + 0.5*erf(erfParams{i}(1)*xProj+erfParams{i}(2));

        %if elecCombAxes(combIDs(i)) == -1 %signifies primary-alone stimulation
        if length(elecCombs{combIDs(i)}) == 1 && pAmps(i) < 0 %cathodal primary-alone stimulation
            for j = 1:6
                axes(outerAxes{j})
                hold on
                if plotLogScale
                    plot(log(xProj), projection, 'k-', 'LineWidth', 2)
                else
                    plot(xProj, projection, 'k-', 'LineWidth', 2)
                end
                hold off
            end
            %threshCurvePAlone = thresholds(i);
            %threshCurvePAloneStd = threshStds(i);
        elseif pAmps(i) < 0 && length(elecCombs{combIDs(i)}) == 2 %cathodal primary, single secondary
            sElecInd = find(sAmps{i});
            axes(outerAxes{secondaryAxes(sElecInd)})
            hold on
            if plotLogScale
                current = plot(log(xProj), projection);
            else
                current = plot(xProj, projection);
            end
            for j = 1:nRelAmps
                if abs((actualRelAmps{i}(sAmps{i}~=0) - relAmps(j))/relAmps(j)) < 0.1; %actual and target relative amplitudes within 10% error
                    set(findobj(current,'Type','line'),'Color', colorScheme(j,:), 'LineWidth', 2)
                    threshCurves{sElecInd}(j) = thresholds(i);
                    threshCurveStds{sElecInd}(j) = threshStds(i);
                    threshCurveActualRelAmps{sElecInd}(j) = actualRelAmps{i}(sAmps{i}~=0);
                end
            end
            hold off
        end
    end
end


xTickLocationsAll = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 2.5 3 3.5 4 4.5 5];
xTickLabelsAll = {[], '0.2', [], '0.4', [], '0.6', [], '0.8', [], '1', [], '2', [], '3', [], '4', [], '5'};

iLabel = 0;
for i = 1:length(xTickLocationsAll)
    if xTickLocationsAll(i) >= curvePlotXLims(1) && xTickLocationsAll(i) <= curvePlotXLims(2)
        iLabel = iLabel + 1;
        xTickLabels{iLabel} = xTickLabelsAll{i}; %#ok<AGROW>
    end
end

xTickLocations = log(xTickLocationsAll(xTickLocationsAll >= curvePlotXLims(1) & xTickLocationsAll <= curvePlotXLims(2)));

for i = 1:6
    if plotLogScale
        %xlabel('log(current amplitude, \muA)')
        
        set(outerAxes{i}, 'xLim', log(curvePlotXLims), 'xTickLabel', xTickLabels, 'xTick', xTickLocations)
        axes(outerAxes{i})
        
        xlabel('current amplitude (\muA)')
    else
        set(outerAxes{i}, 'xLim', curvePlotXLims)
        axes(outerAxes{i})
        xlabel('current amplitude, \muA')
    end
    ylabel('response probability')
end

%% estimates expected value of threshold of primary alone by least-squares linear regression of each
% set of thresholds and then averaging
% also determines "slope" value for effect of each secondary (in terms of log(thresh))

% changed to reflect better linear fit in terms of absolute (vs relative) amplitude on
% secondary 2010-04-07
%
%

estimatedThreshes = zeros(6,1);
threshShiftSlopes = zeros(6, 1); %currently not used in later code
for i = 1:6
    xVals = threshCurves{i}.*threshCurveActualRelAmps{i};
    yVals = threshCurves{i};
    
    for j = nRelAmps:-1:1
        if yVals(j) == inf
            disp('warning: not all expected primary+secondary thresholds were found')
            xVals(j) = [];
            yVals(j) = [];
        end
    end
    
    if length(xVals)>1 %makes sure that at least 2 values exist before fitting line
        [threshShiftSlopes(i) estimatedThreshes(i)] = fitLine(xVals, yVals);
    else
        threshShiftSlopes(i) = inf; estimatedThreshes(i) = inf;
    end
end

for i = 6:-1:1
    if estimatedThreshes(i) == inf
        estimatedThreshes(i) = [];
        disp('warning: not all secondary electrodes had enough data to contribute to estimate of primary-alone threshold')
    end
end
estimatedPrimThresh = mean(estimatedThreshes);


%% plots of thresholds and comparison to linear model


insetAxes = cell(6,1);
for ii = 1:6 %corresonds with sElecs
    outerAxesPos = get(outerAxes{secondaryAxes(ii)}, 'position');
    insetPos = [outerAxesPos(1)+outerAxesPos(3)*0.65 outerAxesPos(2)+outerAxesPos(4)*0.1 outerAxesPos(3)*0.35 outerAxesPos(4)*0.3];
    insetAxes{ii} = axes('position', insetPos);
    
    plotLinearityDataVsModel(thresholds, threshStds, details, ii)
    axis equal
    set(gca, 'ylim', insetYLim)
    %axis equal
end


%% multiple secondaries

tripletPatternsBin = zeros(nPatterns, 1);
pairPatternsBin = zeros(nPatterns, 1);
for i = 1:nPatterns
    if length(elecs{i}) == 3
        tripletPatternsBin(i) = 1;
    elseif length(elecs{i}) == 2
        pairPatternsBin(i) = 1;
    end
end
tripletPatternsInd = find(tripletPatternsBin);
pairPatternsInd = find(pairPatternsBin);
nTriplets = sum(tripletPatternsBin);
nPairs = sum(pairPatternsBin);

if nTriplets == 0
    return
end

secondaryCombs = cell(nTriplets, 1); %pattern numbers of individual primary+2 secondaries combination,
%primary+1 secondary, and primary+other secondary


%determines patterns corresponding to primary + single secondary for each secondary in triplet
excludeTriplet = zeros(nTriplets, 1); %true when corresponding to "excluded" pattern numbers
for i=1:nTriplets
    secondaryCombs{i}(1) = tripletPatternsInd(i); %pattern number of triplet
    sAmpsInd = find(sAmps{tripletPatternsInd(i)});
    for j=1:2
        foundMatch = 0;
        sAmp = sAmps{tripletPatternsInd(i)}(sAmpsInd(j));
        for k=1:nPairs
            if (abs(sAmps{pairPatternsInd(k)}(sAmpsInd(j)) - sAmp) < ampMatchTol) && ~foundMatch
                foundMatch = 1;
                secondaryCombs{i}(j+1) = pairPatternsInd(k);
            elseif (abs(sAmps{pairPatternsInd(k)}(sAmpsInd(j)) - sAmp) < ampMatchTol) && foundMatch
                error('more than one individual secondary pattern matched one of the secondary stim amps in the triplet')
            end
        end
        if foundMatch == 0
            error('no individual secondary pattern matched one of the secondary stim amps in the triplet')
        end
    end
    for j=1:3
        if any(excludePatterns == patternNos(secondaryCombs{i}(j)))
            excludeTriplet(i) = 1;
        end
    end
end


details.secondaryCombsIncluded = cell(nTriplets - sum(excludeTriplet), 1);
includeI = 0;
for i = 1:nTriplets
    if ~excludeTriplet(i)
        includeI = includeI + 1;
        details.secondaryCombsIncluded{includeI} = secondaryCombs{i};
    end
end


return


%% makes slider plot: needs to be updated to reflect log vs. not-log
sliderFig = figure;
slider = make_loop_slider_list(1,1,nTriplets);
plotColors = hsv(3);

while ishandle(sliderFig)
    i = round(get(slider,'Value'));
    cla
    threshShift = zeros(2,1);
    hold on
    for j = 1:3
        xProj = log(curvePlotXLims(1):0.01:curvePlotXLims(2));
        projection = 0.5 + 0.5*erf(erfParams{secondaryCombs{i}(j)}(1)*xProj+erfParams{secondaryCombs{i}(j)}(2));

        plot(allData{secondaryCombs{i}(j)}(1,:), allData{secondaryCombs{i}(j)}(2,:),...
            '.', 'markerFaceColor', plotColors(j,:), 'markerEdgeColor', plotColors(j,:));
        plot(xProj, projection, 'Color', plotColors(j,:), 'LineWidth', 2)
        
        if j == 2 || j == 3
            threshShift(j-1) = thresholds(secondaryCombs{i}(j)) - estimatedPrimThresh;
        end
    end
    
    %plots primary-alone curve
    param_1 = erfParams{primAloneInd}(1); %for now, use slope of actual curve
    param_2 = -estimatedPrimThresh*param_1;
    
    xProj = log(curvePlotXLims(1):0.01:curvePlotXLims(2));
    projection = 0.5 + 0.5*erf(param_1*xProj+param_2);
    current = plot(xProj, projection);
    set(findobj(current,'Type','line'),'Color', [0 0 0], 'LineWidth', 1)
    
    %plots expected combination of secondaries curve
    param_1 = erfParams{secondaryCombs{i}(1)}(1); %slope
    param_2 = -(estimatedPrimThresh+threshShift(1)+threshShift(2))*param_1;
    projection = 0.5 + 0.5*erf(param_1*xProj+param_2);
    plot(xProj, projection, 'Color', plotColors(1,:), 'LineWidth', 1)
    
    
    hold off
    set(gca, 'xLim', [-1.6 -0.8])
    uiwait;
end
clear all

pathToData = '/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data003';

relAmps = [-2 -1 -0.5 0.5 1 2];
neuronID = 63;
pElec = 7;
patternNos = [1:36 253 254 256 258 260 262 264];

[thresholds threshStds estimatedPrimThresh details] =...
    generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID, relAmps,...
    'constrainSlopes', 0, 'plotLogScale', 1, 'recalcAll', 0);

% get secondary-alone thresholds
sAloneThresh = zeros(length(details.sElecs), 1);
sAloneThreshStd = zeros(length(details.sElecs), 1);
for i = 1:length(patternNos)
    if details.pAmps(i) == 0
        for j = 1:length(details.sElecs)
            if details.sAmps{i}(j) ~= 0
                sAloneThresh(j) = thresholds(i);
                sAloneThreshStd(j) = threshStds(i);
            end
        end
    end
end


sumPlot = figure;
hold on
for i = 1:6
    if sAloneThresh(i) ~= inf
        plot(sAloneThresh(i), details.threshShiftSlopes(i), 'k.')
    else
        plot(1.5, details.threshShiftSlopes(i), 'k*')
    end
end
hold off

clear thresholds threshStds estimatedPrimThresh details

%%


pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';
relAmps = [-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1];
neuronID = 800;
pElec = 49;

patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317 319 321 323 325 327 329];
hold off

[thresholds threshStds estimatedPrimThresh details] =...
    generateClusterStimSummaryPlots(pathToData, patternNos, pElec, neuronID, relAmps,...
    'constrainSlopes', 0, 'plotLogScale', 1, 'recalcAll', 0);

% get secondary-alone thresholds
sAloneThresh = zeros(length(details.sElecs), 1);
sAloneThreshStd = zeros(length(details.sElecs), 1);
for i = 1:length(patternNos)
    if details.pAmps(i) == 0
        for j = 1:length(details.sElecs)
            if details.sAmps{i}(j) ~= 0
                sAloneThresh(j) = thresholds(i);
                sAloneThreshStd(j) = threshStds(i);
            end
        end
    end
end

figure
figure(sumPlot)
hold on
for i = 1:6
    if sAloneThresh(i) ~= inf
        plot(sAloneThresh(i), details.threshShiftSlopes(i), 'b.')
    else
        plot(1.5, details.threshShiftSlopes(i), 'b*')
    end
end
hold off
%%
clear thresholds threshStds estimatedPrimThresh details

%%
pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008';

relAmps = [-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1];
neuronID = 559;

pElec = 44;
patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317 319 321 323 325 327 329]; %only patterns with negative primary

[thresholds, threshStds, estPrimThresh details] = generateClusterStimSummaryPlots(pathToData, patternNos, pElec,...
    neuronID, relAmps, 'plotLogScale', 1, 'constrainSlopes', 0,...
    'threshPlotLims', [0.40 0.65], 'curvePlotXLims', [0.25 1.2], 'redoConstrainedFitting', 0, 'recalcAll', 0);


% get secondary-alone thresholds
sAloneThresh = zeros(length(details.sElecs), 1);
sAloneThreshStd = zeros(length(details.sElecs), 1);
for i = 1:length(patternNos)
    if details.pAmps(i) == 0
        for j = 1:length(details.sElecs)
            if details.sAmps{i}(j) ~= 0
                sAloneThresh(j) = thresholds(i);
                sAloneThreshStd(j) = threshStds(i);
            end
        end
    end
end

figure(sumPlot)

hold on
for i = 1:6
    if sAloneThresh(i) ~= inf
        plot(sAloneThresh(i), details.threshShiftSlopes(i), 'r.')
    else
        plot(1.5, details.threshShiftSlopes(i), 'r*')
    end
end
hold off

%clear thresholds threshStds estimatedPrimThresh details



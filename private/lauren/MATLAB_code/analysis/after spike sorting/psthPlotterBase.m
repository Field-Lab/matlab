function histVals = psthPlotterBase(axesH, binEdges, latenciesMs, varargin)

%% note - latency values == 0 are discarded!!!!!

p = inputParser;

p.addRequired('axesH', @ishandle)
p.addRequired('binEdges', @isnumeric) %edges of each bin, as vector of values
p.addRequired('latenciesMs', @isnumeric)

p.addParamValue('lineColor', [0 0 0], @isnumeric)
p.addParamValue('fillHist', false, @islogical)
p.addParamValue('normalize', false, @islogical)
p.addParamValue('normalizeToPeak', false, @islogical)


p.parse(axesH, binEdges, latenciesMs, varargin{:})

lineColor = p.Results.lineColor;
fillHist = p.Results.fillHist;
normalize = p.Results.normalize;
normToPeak = p.Results.normalizeToPeak;



%% plotting PSTH
axes(axesH)

pathCoords = zeros(2, 2*length(binEdges));

nPulses = length(latenciesMs);

%zero values in latencies signify no response
latenciesForHist = latenciesMs(latenciesMs~=0);

if isempty(latenciesForHist)
    warndlg('no histogram plotted because no spikes present!')
    histVals = [binEdges(1:end-1); zeros(size(binEdges(1:end-1)))];
    return
end

binCounts = histc(latenciesForHist', binEdges);
binCounts = binCounts(1:end-1); %last bin only includes latencies == tempMinEndMs;

histVals = [binEdges(1:end-1); binCounts];

axes(axesH)
%hold on
pathCoords(:, 1) = [binEdges(1) 0];
pathCoords(:, 2) = [binEdges(1) binCounts(1)];
for i = 1:length(binEdges)-2
    pathCoords(:, 3 + 2*(i-1)) = [binEdges(i+1) binCounts(i)];
    pathCoords(:, 4 + 2*(i-1)) = [binEdges(i+1) binCounts(i+1)];
end
pathCoords(:, 2*length(binEdges) - 1) = [binEdges(end) binCounts(end)];
pathCoords(:, 2*length(binEdges)) = [binEdges(end) 0];

if normalize
    pathCoords(2,:) = pathCoords(2,:)/nPulses;
end

if normToPeak
    pathCoords(2,:) = pathCoords(2,:)/max(pathCoords(2,:));
end

plot(pathCoords(1,:), pathCoords(2,:), 'Color', lineColor)

if fillHist
    %fill(pathCoords(1,:), pathCoords(2,:), lineColor, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    fill(pathCoords(1,:), pathCoords(2,:), lineColor);
end

%hold off
set(axesH, 'xLim', [0 binEdges(end)])


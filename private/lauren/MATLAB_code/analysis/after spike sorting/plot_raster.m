function plot_raster(axesH, elecResp, movieNo, options)

if exist('options','var')
    if isfield(options, 'markerSize')
        markerSize = options.markerSize;
    else
        markerSize = 5;
    end
    if isfield(options, 'markerColor')
        markerColor = options.markerColor;
    else
        markerColor = [0 0 0];
    end
    if isfield(options, 'sampleRate')
        sampleRate = options.sampleRate;
    else
        sampleRate = 20000;
    end
    if isfield(options, 'labelYAxis')
        labelYAxis = options.labelYAxis;
    else
        labelYAxis = true;
    end
else
    markerSize = 5;
    markerColor = [0 0 0];
    sampleRate = 20000;
    labelYAxis = true;
end


movieNos = elecResp.stimInfo.movieNos;
movieInd = find(movieNos == movieNo);

if ~elecResp.analysis.finalized(movieInd)
    disp(['warning: movie ' num2str(movieNo) ' of pattern ' num2str(elecResp.stimInfo.patternNo) ' is not locked'])
end

latenciesMs = elecResp.analysis.latencies{movieInd}*1000/sampleRate; %in ms

nPulses = length(latenciesMs);
pulseVector = elecResp.stimInfo.pulseVectors{movieInd};


% determining limits of possible spike times
neuronID = [elecResp.cells.main elecResp.cells.active{movieInd}];
goodChannels = elecResp.cells.goodElecs;
centerChannelIndex = find(goodChannels == elecResp.cells.recElec);
nTemplates = length(neuronID);

templates = cell(nTemplates, 1);
templateMinPos = zeros(nTemplates, 1);

for i = 1:nTemplates
    templates{i} = elecResp.cells.allEIs{elecResp.cells.all == neuronID(i)}(goodChannels, :);
    templateMinPos(i) = find(squeeze(templates{i}(centerChannelIndex,:)) ==...
        min(squeeze(templates{i}(centerChannelIndex,:))));
end
tempMinStart = elecResp.analysis.details.tempOffsetWindow{movieInd}(1) + max(templateMinPos);
tempMinEnd = elecResp.analysis.details.tempOffsetWindow{movieInd}(2) + min(templateMinPos);

tempMinEndMs = tempMinEnd*1000/sampleRate;




%% plotting spike times

axes(axesH)
hold on
for k = 1:nPulses
    if latenciesMs(k) ~= 0
        plot(latenciesMs(k), k, '.', 'markerSize', markerSize,...
            'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
    end
end
hold off
set(axesH, 'xLim', [0 tempMinEnd*1000/sampleRate], 'yLim', [0 length(latenciesMs)+1])
xlabel('time (ms)')
if labelYAxis
    ylabel('trial')
end

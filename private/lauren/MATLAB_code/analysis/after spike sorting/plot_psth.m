function histVals = plot_psth(axesH, elecResp, movieNo, varargin)

% a wrapper function that extracts the needed information from an elecResp struct to plot a psth for
% a specified movie number in the specified axes
%
% uses psthPlotterBase to plot the histogram
%
% needs to be tested


p = inputParser;

p.addRequired('axesH', @ishandle)
p.addRequired('movieNumber')
%p.addRequired('movieNumber', @isnumeric)
p.addRequired('elecResp')
%p.addRequired('elecResp', @isstruct);


p.addParamValue('sampleRate', 20000, @isnumeric)
p.addParamValue('binSize', 0.1, @isnumeric)
p.addParamValue('lineColor', [0 0 0], @isnumeric)

p.parse(axesH, elecResp, movieNo, varargin{:})

sampleRate = p.Results.sampleRate;
binSize = p.Results.binSize;
lineColor = p.Results.lineColor;



%%
movieNos = elecResp.stimInfo.movieNos;

latenciesMs = [];
tempOffsetMax = 0;
for ii = 1:length(movieNo)
    if ~elecResp.analysis.finalized(movieNos == movieNo(ii))
        disp(['warning: movie ' num2str(movieNo(ii)) ' of pattern ' num2str(elecResp.stimInfo.patternNo) ' is not locked'])
    end
    if ~isempty(elecResp.cells.active{movieNos == movieNo(ii)})
        warning(['multiple EIs fitted in analysis of ' elecResp.names.rrs_short_name ', ' num2str(movieNo(ii)) ' - check accuracy of latencies'])
    end
    
    latenciesMs = [latenciesMs; elecResp.analysis.latencies{movieNos == movieNo(ii)}*1000/sampleRate]; %#ok<AGROW> %in ms
    
    tempOffsetMax = max(tempOffsetMax, elecResp.analysis.details.tempOffsetWindow{movieNos == movieNo(ii)}(2));
end
    
nPulses = length(latenciesMs);


% determining limits of possible spike times
%neuronID = [elecResp.cells.main elecResp.cells.active{movieNos == movieNo(ii)}];
neuronID = elecResp.cells.main;
goodChannels = elecResp.cells.goodElecs;
centerChannelIndex = find(goodChannels == elecResp.cells.recElec);
%nTemplates = length(neuronID);

template = squeeze(elecResp.cells.allEIs{elecResp.cells.all == neuronID}(elecResp.cells.recElec, :));
templateMinPos = find(template == min(template));
tempMinEnd = tempOffsetMax + templateMinPos;

%templates = cell(nTemplates, 1);
%templateMinPos = zeros(nTemplates, 1);
%for i = 1:nTemplates
%    templates{i} = elecResp.cells.allEIs{elecResp.cells.all == neuronID(i)}(goodChannels, :);
%    templateMinPos(i) = find(squeeze(templates{i}(centerChannelIndex,:)) ==...
%        min(squeeze(templates{i}(centerChannelIndex,:))));
%end
%tempMinEnd = elecResp.analysis.details.tempOffsetWindow{movieInd}(2) + min(templateMinPos);

tempMinEndMs = tempMinEnd*1000/sampleRate;

%final bin should end at or before tempMinEndMs to avoid clipping effects
binEdges = 0:binSize:tempMinEndMs+1;


%% plotting PSTH

histVals = psthPlotterBase(axesH, binEdges, latenciesMs, 'lineColor', lineColor);

xlabel('time (ms)')
ylabel('number of spikes')

set(axesH, 'yLim', [0 nPulses])


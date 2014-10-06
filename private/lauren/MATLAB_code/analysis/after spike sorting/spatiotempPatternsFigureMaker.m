function spatiotempPatternsFigureMaker(pathToData, neuronIDs, patternNos, patternTimes, movieNo, varargin)


%neurons/patterns will be displayed from bottom to top

%neuronIDs         neuron targets should correspond with pattern numbers (matches colors)
%patternNos        should match actual order of pulses, determined by corresponding pattern times!!!
%patternTimes      times at which corresponding patternNos are played (in ms)
%electrodes        if specified, should correspond with patternNos


p = inputParser;

p.addRequired('pathToData', @ischar)
p.addRequired('neuronIDs', @isnumeric)
p.addRequired('patternNos', @isnumeric)
p.addRequired('patternTimes', @isnumeric)
p.addRequired('movieNo', @isnumeric)

colorSchemeDef(1,:) = [27 117 187]/256;
colorSchemeDef(2,:) = [190 30 45]/256;
colorSchemeDef(3,:) = [41 180 115]/256;
colorSchemeDef(4,:) = [0.8 0.5 0];

p.addParamValue('colorScheme', colorSchemeDef, @isnumeric)
p.addParamValue('xLimits', [], @isnumeric)
p.addParamValue('electrodes', [], @isnumeric) % this is required for simultaneous stimulation patterns!!

p.parse(pathToData, neuronIDs, patternNos, patternTimes, movieNo, varargin{:})

colorScheme = p.Results.colorScheme;
xLimits = p.Results.xLimits;
electrodes = p.Results.electrodes;

if isempty(xLimits)
    if length(patternNos) > 1 %not simultaneous stim
        patternTimesInOrder = sort(patternTimes);
        xLimits = [0 (patternTimesInOrder(2)-patternTimesInOrder(1))*length(patternTimes)];
    else
        xLimits = [0 1];
    end
end

prevDir = pwd;

cd(pathToData)


%% extract data
latencies = cell(length(neuronIDs), length(patternNos)); %stores latency information from elecResps
pulses = cell(length(patternNos), 1);
pulseAmps = zeros(length(patternNos), 1);
electrodesAll = cell(length(patternNos), 1);
for j = 1:length(patternNos) % subsequent patterns
    for i = 1:length(neuronIDs) %subsequent neurons
        temp = load(['elecResp_n' num2str(neuronIDs(i)) '_p' num2str(patternNos(j))]);
        elecResp = temp.elecResp;
        if ~elecResp.analysis.finalized(movieNo)
            disp(['warning: movie ' num2str(movieNo) ' of pattern ' num2str(patternNos(i,j)) ' is not locked'])
        end
        movieIndex = find(elecResp.stimInfo.movieNos == movieNo);
        latencies{i,j} = elecResp.analysis.latencies{movieIndex};
    end
    pulses{j} = elecResp.stimInfo.pulseVectors{movieIndex}; % elec x time (1) vs. amp (2) x samples
    pulseAmps(j) = elecResp.stimInfo.stimAmps(movieIndex);
    electrodesAll{j} = elecResp.stimInfo.electrodes;
end


%% check for refractory period violations (refractory period of 1 ms)

for i = 1:length(neuronIDs)
    allLatencies = [];
    for j = 1:length(patternNos)
        allLatencies = [latencies{i,j}/20 + patternTimes(j)*(latencies{i,j}~=0), allLatencies];
    end
    for j = 1:size(allLatencies, 1) %goes through trials
        latenciesThisTrial = allLatencies(j, allLatencies(j,:) ~= 0);
        latencyIntervals = diff(sort(latenciesThisTrial));
        if any(latencyIntervals < 0.5)
            warndlg('There is a refractory period (0.5 ms) violation in this data')
        end
    end
end

%% make figure

figure('position', [100 100 300 300])
    
% plotting pulses

pulseAxes = axes('position', [0.05 0.05 0.2*xLimits(2) 0.2]);
labelAxes = axes('position', [0.2*xLimits(2)+0.05 0.05 0.25 0.2]);
for j = 1:length(patternNos) %different patterns
    if size(pulses{j}, 1) == 1 %single electrode
        axes(pulseAxes); hold on
        pulseVector = reshape(pulses{j}(1, :, :), 2, []);
        plot([0, pulseVector(1,:)/1000 + patternTimes(j), xLimits(2)], [2*j - 1, pulseVector(2,:) + 2*j - 1, 2*j - 1],...
            'Color', colorScheme(j,:));
        hold off
        axes(labelAxes); hold on
        text(0.1, 2*j - 1, num2str(electrodesAll{j})); hold off
        pulseYLims = [0 length(patternNos)*2 + 1];

    else %multiple electrodes
        for k = 1:size(pulses{j}, 1)
            elecID = find(electrodesAll{j} == electrodes(k));
            if isempty(elecID)
                error('no match found for one of electrodes in a pattern')
            end
            axes(pulseAxes); hold on
            pulseVector = reshape(pulses{j}(elecID, :, :), 2, []);
            plot([0, pulseVector(1,:)/1000 + patternTimes, xLimits(2)], [2*elecID - 1, pulseVector(2,:) + 2*elecID - 1, 2*elecID - 1],...
                'Color', colorScheme(elecID,:));
            hold off
            axes(labelAxes); hold on
            text(0.1, 2*elecID - 1, num2str(electrodes(elecID)))
            hold off
        end
        pulseYLims = [0 size(pulses{j},1)*2 + 1];
    end
end

set(pulseAxes, 'xLim', xLimits, 'yLim', pulseYLims)
set(labelAxes, 'yLim', get(pulseAxes, 'yLim'))
axes(pulseAxes); xlabel('time (ms)')
axes(labelAxes); axis off

%plotting spikes
for i = 1:length(neuronIDs)
    if size(pulses{j}, 1) == 1 %single electrode
        successRate = sum(latencies{i,i} ~= 0)/length(latencies{i,i});
    else
        successRate = sum(latencies{i,1} ~= 0)/length(latencies{i,1});
    end
    
    %plots neuron IDs and success rates
    axes('position', [0.2*xLimits(2)+0.05 0.275+(i-1)*.175 0.25 0.15])
    axis off
    text(0.1, 0.4, [num2str(neuronIDs(i)), 10,  num2str(successRate)])
    %text(0.1, 0.4, [num2str(neuronIDs(i)), 10,  num2str(successRate), 10, num2str(pulseAmps(i))])
    set(gca, 'xlim', [0 1])
    
    %plots spikes
    axes('position', [0.05 0.275+(i-1)*.175 0.2*xLimits(2) 0.15]); %entire row
    hold on
    for j = 1:length(patternNos)
        for k = 1:length(latencies{i,j}) %plot each spike from ith neuron, jth pattern (raster)
            if latencies{i,j}(k) ~= 0
                % x-value = ms from start of first pulse in series
                plot(latencies{i,j}(k)/20 + patternTimes(j), k, '.', 'markerSize', 5,...
                    'markerFaceColor', colorScheme(i,:), 'markerEdgeColor', colorScheme(i,:))
            end
        end
    end
    hold off

    set(gca, 'xLim', xLimits, 'yLim', [0 length(latencies{i,j})+1],'box', 'on',...
        'xTickLabel', [], 'yTickLabel', [])

end


cd(prevDir)


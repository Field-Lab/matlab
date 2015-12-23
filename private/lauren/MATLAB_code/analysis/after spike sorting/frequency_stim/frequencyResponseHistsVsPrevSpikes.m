function figH = frequencyResponseHistsVsPrevSpikes(elecRespPath, freq, movieInds, nSequence, varargin)

% generates a figure containing histograms of response probability vs.
% pulse # in sequence, and corresonding response curve
%
% ****note that this isn't really meaningful for frequencies that result in
% no gap in time between pulse sequences (e.g. 20 pulses at 20 Hz, repeated
% each second)
%


p = inputParser;

p.addRequired('elecRespPath', @ischar)
p.addRequired('freq', @isnumeric) %frequency of stimulus in Hz
p.addRequired('movieInds', @isnumeric)
p.addRequired('nSequence', @isnumeric) %number of pulses in a single repetition

p.addParamValue('markerSize', 5, @isnumeric) %dots in response curve
p.addParamValue('panelH', [], @ishandle) %if only one handle is given, only plots 
p.addParamValue('respCurvePanelH', [], @ishandle)


p.parse(elecRespPath, freq, movieInds, nSequence, varargin{:})


markerSize = p.Results.markerSize;
panelH = p.Results.panelH;
respCurvePanelH = p.Results.respCurvePanelH;

%%

nMovies = length(movieInds);

plotColors = hsv(nMovies);

if isempty(panelH)
    figH = figure('position', [100 100 400 300]);
end

for jj = 1:nMovies
    
    movieIndTemp = movieInds(jj);
    
    [latencies successes] = extractFrequencyAnalysis(elecRespPath, movieIndTemp, nSequence);
    
    nReps = size(latencies,2);
    
    
    
    % determine probability of response as a function of number of
    % preceding successes
    trialCount = zeros(1,nSequence); %keeps track of number of times this situation occurs
    succCount = zeros(1,nSequence); %keeps track of number of successes when this occurs
    
    for ii = 1:nReps
        currentRun = 1;
        for kk = 1:nSequence
            trialCount(currentRun) = trialCount(currentRun) + 1;
            if successes(kk,ii)
                succCount(currentRun) = succCount(currentRun) + 1;
                currentRun = currentRun + 1;
            else %reset length of spike sequence to 1
                currentRun = 1;
            end
        end
    end
    
    probLength = find(trialCount==0,1) - 1;
    if isempty(probLength)
        probLength = nSequence;
    end
    
    probSequence = succCount(1:probLength)./trialCount(1:probLength);
    

    
    
    % generate sequence of values for histogram-shaped plot
    xToPlot = [0.5];
    yToPlot = [0];
    for ii = 1:probLength
        xToPlot = [xToPlot ii-0.5 ii+0.5];
        yToPlot = [yToPlot probSequence(ii) probSequence(ii)];
    end
    xToPlot = [xToPlot probLength+0.5];
    yToPlot = [yToPlot 0];
    
    if ~isempty(panelH)
        axes('parent', panelH, 'position', [0.2 1-(jj+0.5)/(nMovies+2) 0.6 0.8/(nMovies+2)])
    else
        axes('position', [0.1 1-(jj+0.5)/(nMovies+2) 0.3 0.8/(nMovies+2)])
    end
    
    hold on
    p = patch(xToPlot, yToPlot, [0.8 0.8 0.8]);
    set(p, 'edgeColor', 'none')
    
    plot(xToPlot, yToPlot, 'color', plotColors(jj,:))
    hold off
    
    %ylabel([num2str(stimAmp) ' µA'])
    set(gca,'ylim', [-0.1 1.1], 'xlim', [0 nSequence + 1], 'ytick', [0 1])
    
    if jj ~= nMovies
        set(gca, 'xticklabel', [])
    else
        xlabel('spike number')
    end
    
    if jj == 1
        title('Resp. Prob. vs spike in sequence')
    end
    
end

%% plot corresponding response curve

if isempty(panelH) || ~isempty(respCurvePanelH)
    
    temp = load(elecRespPath);
    
    shortName = temp.elecResp.names.rrs_short_name;
    fIndex = strfind(shortName, '_f');
    shortName = shortName(1:fIndex-1);
    underIndex = strfind(shortName, '_');
    for ii = 1:length(underIndex)
        shortName = [shortName(1:underIndex(ii)-1) ' ' shortName(underIndex(ii)+1:end)];
    end
    
    if isempty(respCurvePanelH)
        axes('position', [0.5 0.8 0.5 0.15], 'xlim', [0 1], 'ylim', [0 1])
    else
        axes('parent', respCurvePanelH, 'position', [0.05 0.8 0.9 0.15], 'xlim', [0 1], 'ylim', [0 1])
    end
    text(0.01, 0.5, [shortName, 10, 'frequency = ', num2str(freq) ' Hz'])
    axis off
    
    
    data = zeros(2, length(temp.elecResp.stimInfo.movieNos));
    data(1,:) = temp.elecResp.stimInfo.stimAmps;
    data(2,:) = temp.elecResp.analysis.successRates;
    lockedAmps = temp.elecResp.analysis.finalized;
    dotColors = zeros(length(temp.elecResp.stimInfo.movieNos),3);
    dotColors(movieInds,:) = plotColors;
    for i = length(temp.elecResp.stimInfo.stimAmps): -1: 1
        if isempty(temp.elecResp.analysis.type{i})
            data(:,i) = [];
            lockedAmps(i) = [];
            dotColors(i,:) = [];
        end
    end
    data(1,:) = abs(data(1,:));
    
    a = temp.elecResp.analysis.erfParams(1);
    b = temp.elecResp.analysis.erfParams(2);
    
    nAmps = size(data, 2);
    xProj = data(1,1):0.001:data(1,nAmps);
    projection = 0.5 + 0.5*erf(a*xProj+b);
    
    
    if isempty(respCurvePanelH)
        axes('position', [0.55 0.2 0.4 0.5])
    else
        axes('parent', respCurvePanelH, 'position', [0.1 0.2 0.8 0.5])
    end
    
    hold on
    
    title(['threshold = ', num2str(-b/a)])
    for i = 1:size(data, 2)
        if lockedAmps(i)
            plot(data(1,i), data(2,i),'o', 'markerSize', markerSize, 'markerFaceColor', dotColors(i,:), 'markerEdgeColor', [0 0 0]);
        else
            plot(data(1,i), data(2,i),'o', 'markerSize', markerSize, 'markerFaceColor', dotColors(i,:), 'markerEdgeColor', [0.5 0.5 0.5]);
        end
    end
    
    plot(xProj, projection,'k-');
    set(gca, 'ylim', [0 1])
    xlabel('current amplitude (µA)')
    ylabel('response rate')
    hold off

end


clear all

nPulseLimit = 20; %only take into account the first n pulses of pulse sequence for each repetition (only affects pulse sequences that don't fill 1 second movie chunk)

normalizeOver = 1:2; %which of the frequencies to average to determine normalization factor for all frequency thresholds for a given cell

%includes data collected up to (including) 2011-07-14

%% complete data (all curves reach at least 90% or saturate)

for z = 1 %for code-collapsing purposes only!
    
    i = 0;
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n2_p58';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 20; %number of pulses per sequence
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data020/elecResp_n230_p20';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-4/data009/elecResp_n813_p53';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data004/elecResp_n753_p51';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data008/elecResp_n482_p36';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n947_p58';
    elecResps(i).type = 'offPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    elecResps(i).saturateFreq = 160;
    elecResps(i).saturateMovie = 114;
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data011/elecResp_n590_p40';
    elecResps(i).type = 'offPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028/elecResp_n875_p61';
    elecResps(i).type = 'onMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005/elecResp_n616_p45';
    elecResps(i).type = 'onMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data004/elecResp_n186_p10';
    elecResps(i).type = 'offMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    elecResps(i).saturateFreq = [80 160];
    elecResps(i).saturateMovie = [95 102]; %102 could be +/- 1 movie index
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data006/elecResp_n109_p10';
    elecResps(i).type = 'offMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data006/elecResp_n18_p5';
    elecResps(i).type = 'sbc';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data014/elecResp_n797_p47';
    elecResps(i).type = 'sbc';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'complete';
    
    
    %% near complete data/analysis (not all curves reach 90%, but all reach at least 50%)
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-0/data004/elecResp_n437_p30';
    elecResps(i).type = 'onPar';
    elecResps(i).nSequence = 50;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data015/elecResp_n141_p10';
    elecResps(i).type = 'offPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-5/data003/elecResp_n289_p20';
    elecResps(i).type = 'offPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data008/elecResp_n875_p56';
    elecResps(i).type = 'offPar';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data022/elecResp_n274_p16';
    elecResps(i).type = 'onMidg';
    elecResps(i).nSequence = 50;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data005/elecResp_n286_p20';
    elecResps(i).type = 'onMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006/elecResp_n646_p36';
    elecResps(i).type = 'onMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006/elecResp_n648_p44';
    elecResps(i).type = 'offMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data007/elecResp_n621_p42';
    elecResps(i).type = 'offMidg';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2010-10-28-3/data003/elecResp_n242_p15';
    elecResps(i).type = 'sbc';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data004/elecResp_n768_p52';
    elecResps(i).type = 'sbc';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'incomplete';
    
    %% partially complete (not all curves reach 50%, but all reach 30%)
    
    i = i+1;
    elecResps(i).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data007/elecResp_n556_p38';
    elecResps(i).type = 'sbc';
    elecResps(i).nSequence = 20;
    elecResps(i).completeness = 'partial';
    
end

%%
psthParams.markerSize = 5;

plotErfsEachCell = 0;
plotHistEachCellFreq = 160;

frequencies = [5 10 20 40 80 160];

for ii = 1:length(elecResps)
    
    ii
    
    if plotErfsEachCell
        markerSize = 5;
        
        nFreq = length(frequencies);
        figure('position', [100 100 400 300])
        hold on
        
        plotColors = hsv(nFreq);
        legendStrings = cell(1,nFreq);
    end
    
    
    for kk = 1:length(frequencies)
        elecRespPathTemp = [elecResps(ii).path '_f' num2str(frequencies(kk)) '.mat'];
        load(elecRespPathTemp)
        
        %if response curve saturates, only use pre-saturation region for
        %curve-fitting
        if any(elecResps(ii).saturateFreq == frequencies(kk))
            freqIndex = find(elecResps(ii).saturateFreq == frequencies(kk));
            excludeMovies = elecResp.stimInfo.movieNos > elecResps(ii).saturateMovie(freqIndex);
        else
            excludeMovies = false(size(elecResp.stimInfo.movieNos));
        end
        
                
        if nPulseLimit == 0 %don't limit number of pulses used for analysis
            data = zeros(2, length(elecResp.stimInfo.movieNos));
            data(1,:) = elecResp.stimInfo.stimAmps;
            data(2,:) = elecResp.analysis.successRates;
            lockedAmps = elecResp.analysis.finalized;
            movieNos = elecResp.stimInfo.stimAmps;
            for jj = length(elecResp.stimInfo.stimAmps): -1: 1
                if isempty(elecResp.analysis.type{jj}) || excludeMovies(jj)
                    data(:,jj) = [];
                    lockedAmps(jj) = [];
                end
            end
        else %only use first nPulseLimit pulses of each pulse sequence for analysis (only affects pulse sequences with time gaps in between, i.e.
            % sequences that don't fill 1 second movie chunk)
            data = zeros(2,0);
            lockedAmps = [];

            for jj = 1:length(elecResp.stimInfo.movieNos)
                if ~isempty(elecResp.analysis.type{jj}) && ~excludeMovies(jj)
                    data = [data zeros(3,1)]; %#ok<AGROW>
                    lockedAmps = [lockedAmps elecResp.analysis.finalized(jj)]; %#ok<AGROW>
                    data(1,end) = elecResp.stimInfo.stimAmps(jj);

                    if elecResps(ii).nSequence/frequencies(kk) >= 1 %no pauses in pulse sequence
                        successes = elecResp.analysis.latencies{jj} ~= 0;
                        
                        %no pauses in pulse sequence, so just use first nPulseLimit pulses
                        %data(2,end) = sum(successes(1:nPulseLimit))/nPulseLimit;
                        
                        %use all pulses
                        data(2,end) = sum(successes)/length(successes);
                        data(3,end) = length(elecResp.analysis.latencies{jj});
                    else
                        [latencies successes] = extractFrequencyAnalysis(elecRespPathTemp, jj, elecResps(ii).nSequence);
                        data(2,end) = sum(sum(successes(1:nPulseLimit, :)))/(nPulseLimit*size(successes,2));
                        data(3,end) = nPulseLimit*size(successes,2);
                    end
                end
            end
        end
        
        
        
        data(1,:) = abs(data(1,:));
        
        erfParams = erfFitter(data, 2, -1, 'makePlot', 0, 'lockedAmps', lockedAmps);
        %erfParams = erfFitter(data, 2, -1, 'makePlot', 1, 'lockedAmps', lockedAmps);
        
%         
%         dataCompare = zeros(2, length(elecResp.stimInfo.movieNos));
%         dataCompare(1,:) = elecResp.stimInfo.stimAmps;
%         dataCompare(2,:) = elecResp.analysis.successRates;
%         dataCompare(3,:) = length(elecResp.analysis.latencies);
%         for jj = length(elecResp.stimInfo.stimAmps): -1: 1
%             if isempty(elecResp.analysis.type{jj})
%                 dataCompare(:,jj) = [];
%             end
%         end
%         dataCompare(1,:) = abs(dataCompare(1,:));
% 
%         
%         erfFitter(dataCompare, 2, -1, 'makePlot', 1, 'lockedAmps', lockedAmps);
%         xProj = dataCompare(1,1):0.001:dataCompare(1,end);
%         projection =  0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));
%         hold on
%         plot(xProj, projection)
% 
%         keyboard

        
        threshold = -erfParams(2)/erfParams(1);
        
        
%         %determine 50% crossing (direct interpolation between
%         %pair of flanking response probabilities or between 4 response
%         %probablities if there is a back-crossing
        cross50 = [];
        for jj = 1:size(data,2)-1
            if data(2,jj) < 0.5 && data(2,jj+1) >= 0.5
                %linear interpolation
                crosspoint = (0.5 - data(2,jj))/(data(2,jj+1)-data(2,jj))*(data(1,jj+1) - data(1,jj)) + data(1,jj);
                cross50 = [cross50 crosspoint];
            end
        end
        
        if length(cross50) ~= 1
            %figure
            %plot(data(1,:), data(2,:), 'k-')
            warning(['more or less than one 50% crossing in ' elecResps(ii).path '_f' num2str(frequencies(kk)) '.mat'])
            %keyboard
            cross50 = inf;
        end

        elecResps(ii).cross50(kk) = cross50;
        elecResps(ii).threshold(kk) = elecResp.analysis.threshold;
        
        
        
        
        
        if plotErfsEachCell
            xProj = data(1,1):0.001:data(1,end);
            projection = 0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));
            
            hold on
            
            for jj = 1:size(data, 2)
                if lockedAmps(jj)
                    hTemp = plot(data(1,jj), data(2,jj),'.', 'markerSize', markerSize, 'markerEdgeColor', plotColors(kk,:));
                else
                    hTemp = plot(data(1,jj), data(2,jj),'.', 'markerSize', markerSize, 'markerEdgeColor', [0.6 0.6 0.6]);
                end
                set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
            end
            
            hTemp = plot(data(1,:), data(2,:), '--', 'color', plotColors(kk,:));
            set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
            
            hTemp = plot([threshold threshold], [0 0.5], 'color', plotColors(kk,:));
            set(get(get(hTemp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
            
            
            plot(xProj, projection,'-', 'color', plotColors(kk,:));
            legendStrings{kk} = [num2str(frequencies(kk)) 'Hz'];
            if kk == length(frequencies)
                shortName = elecResp.names.rrs_short_name;
                fIndex = strfind(shortName, '_f');
                shortName = shortName(1:fIndex-1);
                underIndex = strfind(shortName, '_');
                for jj = 1:length(underIndex)
                    shortName = [shortName(1:underIndex(jj)-1) ' ' shortName(underIndex(jj)+1:end)];
                end
                
                
                
                legend(legendStrings, 'location', 'SouthEast')
                
                set(gca, 'ylim', [0 1])
                xlabel('current amplitude (µA)')
                ylabel('response rate')
                title(['frequency-dependent responses ', 10, shortName, 10, elecResps(ii).type])
                hold off
            end
        end
    end
end

%% histogram vs. pulse number plots

for ii = 1:length(elecResps)
    for jj = 1:length(plotHistEachCellFreq)

        elecRespPathTemp = [elecResps(ii).path '_f' num2str(plotHistEachCellFreq(jj)) '.mat'];
        load(elecRespPathTemp)
        
        % plot pulse-by-pulse response probabilities (over repetitions
        % of pulse sequence) for a range of stimulus amplitudes
        movieInds = find((elecResp.analysis.successRates < 0.95 & elecResp.analysis.finalized), 6, 'last');
        
        figH = figure('position', [100 100 600 300]);
        panel1 = uipanel('units', 'normalized', 'position', [0 0 0.33 1], 'backgroundColor', [1 1 1], 'borderType', 'none');
        panel2 = uipanel('units', 'normalized', 'position', [0.33 0 0.33 1], 'backgroundColor', [1 1 1], 'borderType', 'none');
        panel3 = uipanel('units', 'normalized', 'position', [0.66 0 0.33 1], 'backgroundColor', [1 1 1], 'borderType', 'none');

        
        %plot response probability vs. pulse number
        frequencyResponseHists(elecRespPathTemp, plotHistEachCellFreq(jj), movieInds, elecResps(ii).nSequence, psthParams, 'panelH', panel1)
        
        % plot response probability vs. number of previous successes
        frequencyResponseHistsVsPrevSpikes(elecRespPathTemp, plotHistEachCellFreq(jj), movieInds, elecResps(ii).nSequence, psthParams,...
            'panelH', panel2, 'respCurvePanelH', panel3)
        
    end
end




%% figure

plotColors = [0.5 0.5 0.5
    0 0 0    
    1 0.3 0.3;
    0.7 0 0;
    0 0 0.8];



figure('position', [100 100 400 600])
threshAxesH = axes('position', [0.1 0.55 0.8 0.4]);
title('normalized thresholds')

cross50AxesH = axes('position', [0.1 0.05 0.8 0.4]);
title('normalized 0.5 crossing')



for ii = 1:length(elecResps)
    %if strcmpi(elecResps(ii).completeness, 'complete')
    %if elecResps(ii).nSequence == 20    
    
        switch lower(elecResps(ii).type)
            case 'onpar'
                plotColor = plotColors(1,:);
            case 'offpar'
                plotColor = plotColors(2,:);
            case 'onmidg'
                plotColor = plotColors(3,:);
            case 'offmidg'
                plotColor = plotColors(4,:);
            case 'sbc'
                plotColor = plotColors(5,:);
            otherwise
                error(['invalid cell type for' elecResps(ii).path ': ' elecResps(ii).type])
        end
        
        %threshold
        axes(threshAxesH)
        hold on
        normThresh = elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver));
        
        plot(normThresh, 'color', plotColor)
        plot(normThresh, 'o', 'markerEdgeColor', plotColor, 'markerFaceColor', plotColor)
        hold off
        
        %50% crossing
        axes(cross50AxesH)
        hold on
        normCross50 = elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver));
        
        xVals = 1:length(frequencies);
        
        plot(xVals(normCross50~=inf), normCross50(normCross50~=inf), 'color', plotColor)
        plot(xVals(normCross50~=inf), normCross50(normCross50~=inf), 'o', 'markerEdgeColor', plotColor, 'markerFaceColor', plotColor)

        
        hold off
        
    %end
    %end
    
end
hold off



%% figure with means +/- SD


%tallies up normalized thresholds
onParThreshs = [];
onPar50Crossings = [];
offParThreshs = [];
offPar50Crossings = [];
onMidgThreshs = [];
onMidg50Crossings = [];
offMidgThreshs = [];
offMidg50Crossings = [];
sbcThreshs = [];
sbc50Crossings = [];

for ii = 1:length(elecResps)
    %if strcmpi(elecResps(ii).completeness, 'complete')
    %if elecResps(ii).nSequence == 20
    switch lower(elecResps(ii).type)
        case 'onpar'
            onParThreshs = [onParThreshs; elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver))];
            onPar50Crossings = [onPar50Crossings; elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver))];
        case 'offpar'
            offParThreshs = [offParThreshs; elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver))];
            offPar50Crossings = [offPar50Crossings; elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver))];
        case 'onmidg'
            onMidgThreshs = [onMidgThreshs; elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver))];
            onMidg50Crossings = [onMidg50Crossings; elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver))];
        case 'offmidg'
            offMidgThreshs = [offMidgThreshs; elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver))];
            offMidg50Crossings = [offMidg50Crossings; elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver))];
        case 'sbc'
            sbcThreshs = [sbcThreshs; elecResps(ii).threshold/mean(elecResps(ii).threshold(normalizeOver))];
            sbc50Crossings = [sbc50Crossings; elecResps(ii).cross50/mean(elecResps(ii).cross50(normalizeOver))];
        otherwise
            error(['invalid cell type for' elecResps(ii).path ': ' elecResps(ii).type])
    end
    
    %end
    %end
end


xVals = 1:length(frequencies);
for ii = 1:size(onParThreshs, 2)
    onParThreshMeans(ii) = mean(onParThreshs(:,ii));
    onPar50CrossingsMeans(ii) = mean(onPar50Crossings(onPar50Crossings(:,ii)~=inf,ii));
    onParThreshSDs(ii) = std(onParThreshs(:,ii));
    onPar50CrossingsSDs(ii) = std(onPar50Crossings(onPar50Crossings(:,ii)~=inf,ii));
    
    offParThreshMeans(ii) = mean(offParThreshs(:,ii));
    offPar50CrossingsMeans(ii) = mean(offPar50Crossings(offPar50Crossings(:,ii)~=inf,ii));
    offParThreshSDs(ii) = std(offParThreshs(:,ii));
    offPar50CrossingsSDs(ii) = std(offPar50Crossings(offPar50Crossings(:,ii)~=inf,ii));
    
    onMidgThreshMeans(ii) = mean(onMidgThreshs(:,ii));
    onMidg50CrossingsMeans(ii) = mean(onMidg50Crossings(onMidg50Crossings(:,ii)~=inf,ii));
    onMidgThreshSDs(ii) = std(onMidgThreshs(:,ii));
    onMidg50CrossingsSDs(ii) = std(onMidg50Crossings(onMidg50Crossings(:,ii)~=inf,ii));
    
    offMidgThreshMeans(ii) = mean(offMidgThreshs(:,ii));
    offMidg50CrossingsMeans(ii) = mean(offMidg50Crossings(offMidg50Crossings(:,ii)~=inf,ii));
    offMidgThreshSDs(ii) = std(offMidgThreshs(:,ii));
    offMidg50CrossingsSDs(ii) = std(offMidg50Crossings(offMidg50Crossings(:,ii)~=inf,ii));
    
    sbcThreshMeans(ii) = mean(sbcThreshs(:,ii));
    sbc50CrossingsMeans(ii) = mean(sbc50Crossings(sbc50Crossings(:,ii)~=inf,ii));
    sbcThreshSDs(ii) = std(sbcThreshs(:,ii));
    sbc50CrossingsSDs(ii) = std(sbc50Crossings(sbc50Crossings(:,ii)~=inf,ii));
end
      



figure('position', [100 100 400 600])
threshAxesH = axes('position', [0.1 0.55 0.8 0.4]);
title('normalized thresholds')

%threshold
hold on
for ii = 1:length(frequencies)
    plot(ii, onParThreshMeans(ii), 'o', 'MarkerEdgeColor', plotColors(1,:))
    plot([ii ii], [onParThreshMeans(ii)+onParThreshSDs(ii) onParThreshMeans(ii)-onParThreshSDs(ii)], 'color', plotColors(1,:))

    plot(ii+0.1, offParThreshMeans(ii), 'o', 'MarkerEdgeColor', plotColors(2,:))
    plot([ii+0.1 ii+0.1], [offParThreshMeans(ii)+offParThreshSDs(ii) offParThreshMeans(ii)-offParThreshSDs(ii)], 'color', plotColors(2,:))
    
    plot(ii+0.2, onMidgThreshMeans(ii), 'o', 'MarkerEdgeColor', plotColors(3,:))
    plot([ii+0.2 ii+0.2], [onMidgThreshMeans(ii)+onMidgThreshSDs(ii) onMidgThreshMeans(ii)-onMidgThreshSDs(ii)], 'color', plotColors(3,:))
    
    plot(ii+0.3, offMidgThreshMeans(ii), 'o', 'MarkerEdgeColor', plotColors(4,:))
    plot([ii+0.3 ii+0.3], [offMidgThreshMeans(ii)+offMidgThreshSDs(ii) offMidgThreshMeans(ii)-offMidgThreshSDs(ii)], 'color', plotColors(4,:))
    
    plot(ii+0.4, sbcThreshMeans(ii), 'o', 'MarkerEdgeColor', plotColors(5,:))
    plot([ii+0.4 ii+0.4], [sbcThreshMeans(ii)+sbcThreshSDs(ii) sbcThreshMeans(ii)-sbcThreshSDs(ii)], 'color', plotColors(5,:))
end


text(1.5, 1.6, ['ON Parasol (n=' num2str(size(onParThreshs,1)) ')'], 'color', plotColors(1,:))
text(1.5, 1.55, ['OFF Parasol (n=' num2str(size(offParThreshs,1)) ')'], 'color', plotColors(2,:))
text(1.5, 1.5, ['ON Midget (n=' num2str(size(onMidgThreshs,1)) ')'], 'color', plotColors(3,:))
text(1.5, 1.45, ['OFF Midget (n=' num2str(size(offMidgThreshs,1)) ')'], 'color', plotColors(4,:))
text(1.5, 1.4, ['SBC (n=' num2str(size(sbcThreshs,1)) ')'], 'color', plotColors(5,:))



hold off

%50% crossing
cross50AxesH = axes('position', [0.1 0.05 0.8 0.4]);
title('normalized 0.5 crossing')

hold on
for ii = 1:length(frequencies)
    plot(ii, onPar50CrossingsMeans(ii), 'o', 'MarkerEdgeColor', plotColors(1,:))
    plot([ii ii], [onPar50CrossingsMeans(ii)+onPar50CrossingsSDs(ii) onPar50CrossingsMeans(ii)-onPar50CrossingsSDs(ii)], 'color', plotColors(1,:))

    plot(ii+0.1, offPar50CrossingsMeans(ii), 'o', 'MarkerEdgeColor', plotColors(2,:))
    plot([ii+0.1 ii+0.1], [offPar50CrossingsMeans(ii)+offPar50CrossingsSDs(ii) offPar50CrossingsMeans(ii)-offPar50CrossingsSDs(ii)], 'color', plotColors(2,:))
    
    plot(ii+0.2, onMidg50CrossingsMeans(ii), 'o', 'MarkerEdgeColor', plotColors(3,:))
    plot([ii+0.2 ii+0.2], [onMidg50CrossingsMeans(ii)+onMidg50CrossingsSDs(ii) onMidg50CrossingsMeans(ii)-onMidg50CrossingsSDs(ii)], 'color', plotColors(3,:))
    
    plot(ii+0.3, offMidg50CrossingsMeans(ii), 'o', 'MarkerEdgeColor', plotColors(4,:))
    plot([ii+0.3 ii+0.3], [offMidg50CrossingsMeans(ii)+offMidg50CrossingsSDs(ii) offMidg50CrossingsMeans(ii)-offMidg50CrossingsSDs(ii)], 'color', plotColors(4,:))
    
    plot(ii+0.4, sbc50CrossingsMeans(ii), 'o', 'MarkerEdgeColor', plotColors(5,:))
    plot([ii+0.4 ii+0.4], [sbc50CrossingsMeans(ii)+sbc50CrossingsSDs(ii) sbc50CrossingsMeans(ii)-sbc50CrossingsSDs(ii)], 'color', plotColors(5,:))
end
hold off















clear all

markerSize = 5;
markerColor = [0 0 0];

%datarun.names.rrs_neurons_path = '/Analysis/lauren/2011-01-11-1/data005-from-data000/data005-from-data000.neurons';
datarun.names.rrs_neurons_path = '/Analysis/lauren/2010-03-05-2/data023-from-data018/data023-from-data018.neurons';
datarun = load_neurons(datarun);

firstTrig = 304;


triggers = datarun.triggers(firstTrig:4:length(datarun.triggers)); %period of flash sequence = 4 triggers (sequence from 1: WGBG)
offset = 0.4; %phase relative to trigger at which to plot (should be positive and < 0.5)

period = max(diff(triggers));

% for 2011-01-11-1/data005-from-data000
%trials 1-25
%chosen.onMidg  = [1787 2297]; %6290
%chosen.offMidg = [7308 7726]; %7726
%chosen.sbc = [5885 4625];
%chosen.onPar   = [5311 7281];
%chosen.offPar  = [317 856]; %856

%nTrials = 25; %just use first 25 trials

chosen.onPar = [3391 3935]; %2971 3391 3935 5371 6935
chosen.offPar = [1472 2012]; %1472 2012 2881 3212
chosen.onMidg = [3076 3241]; %454 2509 2791 3033 3076 3241 3737
chosen.offMidg = [931 4666]; %931 4666 6676 6722 7261
chosen.sbc = [1623 3454];%1623 2793 3454 4624

chosen.onMidg = [454 2509];


%nTrials = length(triggers) - 1;
nTrials = 50;


onParSpikes   = cell(2,nTrials);
offParSpikes  = cell(2,nTrials);
onMidgSpikes  = cell(2,nTrials);
offMidgSpikes = cell(2,nTrials);
sbcSpikes     = cell(2,nTrials);

regionTimes = [0.4913 0.9993 1.4906 1.9986 2.4899] - offset;

if mod(firstTrig,4) == 0 %need to fill out for rest of cases
    regionColors = [0 1 0 -1 0]; %0 = grey, 1 = white, -1 = black
end


for ii = 1:2
    for jj = 1:nTrials
        startTrial = triggers(jj)+offset;
        endTrial = triggers(jj+1)+offset;
        
        %ON parasol
        spikesTemp = datarun.spikes{datarun.cell_ids == chosen.onPar(ii)};
        onParSpikes{ii,jj} = spikesTemp(spikesTemp > startTrial & spikesTemp <= endTrial) - startTrial;
        
        %OFF parasol
        spikesTemp = datarun.spikes{datarun.cell_ids == chosen.offPar(ii)};
        offParSpikes{ii,jj} = spikesTemp(spikesTemp > startTrial & spikesTemp <= endTrial) - startTrial;
        
        %ON midget        
        spikesTemp = datarun.spikes{datarun.cell_ids == chosen.onMidg(ii)};
        onMidgSpikes{ii,jj} = spikesTemp(spikesTemp > startTrial & spikesTemp <= endTrial) - startTrial;
        
        %OFF midget
        spikesTemp = datarun.spikes{datarun.cell_ids == chosen.offMidg(ii)};
        offMidgSpikes{ii,jj} = spikesTemp(spikesTemp > startTrial & spikesTemp <= endTrial) - startTrial;
        
        %SBC
        spikesTemp = datarun.spikes{datarun.cell_ids == chosen.sbc(ii)};
        sbcSpikes{ii,jj} = spikesTemp(spikesTemp > startTrial & spikesTemp <= endTrial) - startTrial;        
    end
end



%% plotting spike times
% ON parasol

figure('position', [300 200 600 200])
aTemp{1} = axes('position', [0.2 0.25 0.7 0.25]);
aTemp{2} = axes('position', [0.2 0.6 0.7 0.25]);
for ii = 1:2
    axes(aTemp{ii})
    hold on
    
    if regionColors(1) == 0
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
    elseif regionColors(1) == -1
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5])
    end
    for jj = 2:5
        if regionColors(jj) == 0
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
        elseif regionColors(jj) == -1
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5], 'edgeColor', 'none')
        end
    end
    
    for jj = 1:nTrials %trials
        for kk = 1:length(onParSpikes{ii,jj}) %spikes
            plot(onParSpikes{ii,jj}(kk), nTrials - jj + 1, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        end
    end

    hold off
    set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1], 'yTick', [1 nTrials], 'yticklabel', [nTrials 1], 'xtick', [0 0.5 1 1.5 2])

    
    %set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1])

    if ii == 2
        title('ON parasol')
        set(gca, 'xticklabel', [])
    else
        xlabel('time (s)')
    end
    ylabel('trial')
    %title(['off parasol ' num2str(chosen.onPar(ii))])
end

%% OFF parasol
figure('position', [300 200 600 200])
aTemp{1} = axes('position', [0.2 0.25 0.7 0.25]);
aTemp{2} = axes('position', [0.2 0.6 0.7 0.25]);
for ii = 1:2
    axes(aTemp{ii})
    hold on
    
    if regionColors(1) == 0
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
    elseif regionColors(1) == -1
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5])
    end
    for jj = 2:5
        if regionColors(jj) == 0
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
        elseif regionColors(jj) == -1
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5], 'edgeColor', 'none')
        end
    end
    
    for jj = 1:nTrials %trials
        for kk = 1:length(offParSpikes{ii,jj}) %spikes
            plot(offParSpikes{ii,jj}(kk), nTrials - jj + 1, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        end
    end

    hold off
    set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1], 'yTick', [1 nTrials], 'yticklabel', [nTrials 1], 'xtick', [0 0.5 1 1.5 2])
    if ii == 2
        set(gca, 'xticklabel', [])
        title('OFF parasol')
    else
        xlabel('time (s)')
    end
    ylabel('trial')
    %title(['off parasol ' num2str(chosen.offPar(ii))])
end

%% ON midget

figure('position', [300 200 600 200])
aTemp{1} = axes('position', [0.2 0.25 0.7 0.25]);
aTemp{2} = axes('position', [0.2 0.6 0.7 0.25]);
for ii = 1:2
    axes(aTemp{ii})
    hold on
    
    if regionColors(1) == 0
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
    elseif regionColors(1) == -1
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5])
    end
    for jj = 2:5
        if regionColors(jj) == 0
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
        elseif regionColors(jj) == -1
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5], 'edgeColor', 'none')
        end
    end
    
    for jj = 1:nTrials %trials
        for kk = 1:length(onMidgSpikes{ii,jj}) %spikes
            plot(onMidgSpikes{ii,jj}(kk), nTrials - jj + 1, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        end
    end

    hold off
    set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1], 'yTick', [1 nTrials], 'yticklabel', [nTrials 1], 'xtick', [0 0.5 1 1.5 2])
    if ii == 2
        set(gca, 'xticklabel', [])
        title('ON midget')
    else
        xlabel('time (s)')
    end
    ylabel('trial')
    %title(['off parasol ' num2str(chosen.onMidg(ii))])
end


%% OFF midget

figure('position', [300 200 600 200])
aTemp{1} = axes('position', [0.2 0.25 0.7 0.25]);
aTemp{2} = axes('position', [0.2 0.6 0.7 0.25]);
for ii = 1:2
    axes(aTemp{ii})
    hold on
    
    if regionColors(1) == 0
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
    elseif regionColors(1) == -1
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5])
    end
    for jj = 2:5
        if regionColors(jj) == 0
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
        elseif regionColors(jj) == -1
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5], 'edgeColor', 'none')
        end
    end
    
    for jj = 1:nTrials %trials
        for kk = 1:length(offMidgSpikes{ii,jj}) %spikes
            plot(offMidgSpikes{ii,jj}(kk), nTrials - jj + 1, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        end
    end

    hold off
    set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1], 'yTick', [1 nTrials], 'yticklabel', [nTrials 1], 'xtick', [0 0.5 1 1.5 2])
    if ii == 2
        set(gca, 'xticklabel', [])
        title('OFF midget')
    else
        xlabel('time (s)')
    end
    ylabel('trial')
    %title(['off parasol ' num2str(chosen.offMidg(ii))])
end

%% SBC

figure('position', [300 200 600 200])
aTemp{1} = axes('position', [0.2 0.25 0.7 0.25]);
aTemp{2} = axes('position', [0.2 0.6 0.7 0.25]);
for ii = 1:2
    axes(aTemp{ii})
    hold on
    
    if regionColors(1) == 0
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
    elseif regionColors(1) == -1
        patch([0 regionTimes(1) regionTimes(1) 0], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5])
    end
    for jj = 2:5
        if regionColors(jj) == 0
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.75 0.75 0.75], 'edgeColor', 'none')
        elseif regionColors(jj) == -1
            patch([regionTimes(jj-1) regionTimes(jj) regionTimes(jj) regionTimes(jj-1)], [0 0 nTrials+1 nTrials+1], [0.5 0.5 0.5], 'edgeColor', 'none')
        end
    end
    
    for jj = 1:nTrials %trials
        for kk = 1:length(sbcSpikes{ii,jj}) %spikes
            plot(sbcSpikes{ii,jj}(kk), nTrials - jj + 1, '.', 'markerSize', markerSize,...
                'markerFaceColor', markerColor, 'markerEdgeColor', markerColor)
        end
    end

    hold off
    set(gca, 'xLim', [0 period], 'yLim', [0 nTrials+1], 'yTick', [1 nTrials], 'yticklabel', [nTrials 1], 'xtick', [0 0.5 1 1.5 2])
    if ii == 2
        set(gca, 'xticklabel', [])
        title('SBC')
    else
        xlabel('time (s)')
    end
    ylabel('trial')
    %title(['off parasol ' num2str(chosen.sbc(ii))])
end


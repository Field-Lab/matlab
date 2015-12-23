function plotSpatioTempProbeResults_1elec(prePulse, preElecInd, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)



markerSize = 4;


offsets(offsets == 0) = [];
offsetsMs = offsets/20;


nOffsets = length(offsets);

plotColors(1,:) = [90 156 0]/255; %pale grass
plotColors(2,:) = [255 124 59]/255; %salmon
plotColors(3,:) = [101 52 255]/255; %purple
plotColors(4,:) = [52 198 247]/255; %aqua
plotColors(5,:) = [238 55 128]/255; %calm magenta

%load target cell ei from one of the elecResps
temp = load([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '.mat']);
threshIso = temp.elecResp.analysis.threshold;
threshStdIso = temp.elecResp.analysis.threshStd;

for ii = 1:64
    if ~(ii==9||ii==25||ii==57)
        eiData(ii) = max(abs(temp.elecResp.cells.mainEI(ii,:))); %#ok<AGROW>
    end
end
clear temp



ii = preElecInd;
figure('position', [100 100 620 300], 'color', [1 1 1])

%ei plot
eiAxes = axes('position', [0.03 0.43 0.26 0.53]);
plotEi61('', neuronID, 'axesH', eiAxes, 'eiData', {eiData}, 'markElecs', [stimElec prePulse(ii).elec], 'manualScale', 0.005)

%thresholds plot
threshAxes = axes('position', [0.48 0.2 0.48 0.73]);
nAmps = length(prePulse(ii).amps);

hold on

%plot threshold for isolated stim pulse
plot([0 offsetsMs(end)+5.5], [threshIso threshIso], 'color', [0.8 0.8 0.8])

plot(offsetsMs(end)+5, threshIso, 'o', 'markerFaceColor', 'none', 'markerEdgeColor', [0 0 0], 'markerSize', markerSize)
plot([offsetsMs(end)+5 offsetsMs(end)+5], [threshStdIso, -1*threshStdIso] + threshIso,...
    '-', 'color', [0 0 0])

for kk = 1:nAmps
    
    plotColorSucc = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(4,:);
    plotColorFail = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(5,:);
    
    
    %plot lines connecting data points
    psThresh = prePulse(ii).threshPostSucc(:,kk);
    pfThresh = prePulse(ii).threshPostFail(:,kk);
    
    %thresholds outside measured amp range are plotted at limits of amp range
    %psThresh(isinf(psThresh) & psThresh > 0) = ampRange(2);
    %psThresh(isinf(psThresh) & psThresh < 0) = ampRange(1);
    %pfThresh(isinf(pfThresh) & pfThresh > 0) = ampRange(2);
    %pfThresh(isinf(pfThresh) & pfThresh < 0) = ampRange(1);
    
    if all(isfinite(psThresh))
        plot(offsetsMs, psThresh, 'color', plotColorSucc)
    else
        for nn = 1:length(offsetsMs)-1
            if isfinite(psThresh(nn)) && isfinite(psThresh(nn)) %current and next points are valid
                plot(offsetsMs(nn:nn+1), psThresh(nn:nn+1), 'color', plotColorSucc)
            end
        end
    end
    if all(isfinite(pfThresh)) && prePulse(ii).analyzedPreResponse(kk) == 1
        plot(offsetsMs, pfThresh, 'color', plotColorFail)
    elseif ~all(isfinite(pfThresh))
        for nn = 1:length(offsetsMs)-1
            if isfinite(pfThresh(nn)) && isfinite(pfThresh(nn)) %current and next points are valid
                plot(offsetsMs(nn:nn+1), pfThresh(nn:nn+1), 'color', plotColorFail)
            end
        end
    end
    
    for jj = 1:nOffsets
        plottedPostSucc = false;
        plottedPostFail = false;
        
        %plots post-success thresholds
        if prePulse(ii).threshPostSucc(jj,kk)~=0 && isfinite(prePulse(ii).threshPostSucc(jj,kk))
            plot(offsetsMs(jj), prePulse(ii).threshPostSucc(jj,kk), 'o',...
                'markerEdgeColor', plotColorSucc, 'markerFaceColor', 'none', 'markerSize', markerSize)
            plot([offsetsMs(jj) offsetsMs(jj)],...
                [prePulse(ii).threshStdPostSucc(jj,kk), -1*prePulse(ii).threshStdPostSucc(jj,kk)] + prePulse(ii).threshPostSucc(jj,kk),...
                'color', plotColorSucc)
            plottedPostSucc = true;
        elseif isinf(prePulse(ii).threshPostSucc(jj,kk)) %plots cases where threshold is outside measured amplitudes
            if prePulse(ii).threshPostSucc(jj,kk) > 0
                plot(offsetsMs(jj), ampRange(2), '^', 'markerEdgeColor', plotColorSucc, 'markerSize', markerSize)
            elseif prePulse(ii).threshPostSucc(jj,kk) < 0
                plot(offsetsMs(jj), ampRange(1), 'v', 'markerEdgeColor', plotColorSucc, 'markerSize', markerSize)
            end
            plottedPostSucc = true;
        elseif prePulse(ii).threshPostSucc(jj,kk)==0
            keyboard
        end
        
        if prePulse(ii).analyzedPreResponse(kk) == 1
            if prePulse(ii).threshPostFail(jj,kk)~=0 && isfinite(prePulse(ii).threshPostFail(jj,kk))
                plot(offsetsMs(jj), prePulse(ii).threshPostFail(jj,kk), 'o',...
                    'markerEdgeColor', plotColorFail, 'markerFaceColor', 'none', 'markerSize', markerSize)
                plot([offsetsMs(jj) offsetsMs(jj)],...
                    [prePulse(ii).threshStdPostFail(jj,kk), -1*prePulse(ii).threshStdPostFail(jj,kk)] + prePulse(ii).threshPostFail(jj,kk),...
                    'color', plotColorFail)
                plottedPostFail = true;
            elseif isinf(prePulse(ii).threshPostFail(jj,kk)) %plots cases where threshold is outside measured amplitudes
                if prePulse(ii).threshPostFail(jj,kk) > 0
                    plot(offsetsMs(jj), ampRange(2), '^', 'markerEdgeColor', plotColorFail, 'markerSize', markerSize)
                elseif prePulse(ii).threshPostFail(jj,kk) < 0
                    plot(offsetsMs(jj), ampRange(1), 'v', 'markerEdgeColor', plotColorFail, 'markerSize', markerSize)
                end
                plottedPostFail = true;
            elseif prePulse(ii).threshPostFail(jj,kk)==0
                keyboard
            end
        end
        
        if ~plottedPostFail && ~prePulse(ii).analyzedPreResponse(kk) == 1
            disp(['no data point plotted for elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' num2str(prePulse(ii).amps(kk)) ' (post failure)'])
        end
        if ~plottedPostSucc
            disp(['no data point plotted for elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' num2str(prePulse(ii).amps(kk)) ' (post success)'])
        end
    end
end

set(threshAxes, 'xlim', [0 offsetsMs(end)+5.5], 'ylim', threshYLims, 'xtick', 0:2:offsetsMs(end))
hold off
ylabel('threshold (µA)')
xlabel('delay between prepulse and stim pulse (ms)')


%legends
legendAxes = axes('position', [0.015 0.03 0.275 0.37]);
hold on
for kk = 1:nAmps
    plotColorSucc = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(4,:);
    plotColorFail = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(5,:);
    
    
    if prePulse(ii).analyzedPreResponse(kk) == 1
        plot(0.7-0.2*kk, 0.5, 's', 'MarkerSize', 5,...
            'MarkerEdgeColor', plotColorFail, 'MarkerFaceColor', plotColorFail)
    end
    plot(0.7-0.2*kk, 0.4, 's', 'MarkerSize', 5,...
        'MarkerEdgeColor', plotColorSucc, 'MarkerFaceColor', plotColorSucc)
    
    text(0.7-0.2*kk, 0.6, [num2str(prePulse(ii).amps(kk)) ' µA'], 'rotation', 90)
    
    %mean response probability to prepulse
    text(0.78-0.2*kk, 0.6, ['(' num2str(0.01*round(100*prePulse(ii).meanPreRespProb(kk))) ')'], 'rotation', 90)
end
text(0.7, 0.5, 'after no response')
text(0.7, 0.4, 'after response')

set(legendAxes, 'ylim', [0.2 1], 'xlim', [0 1])
hold off
axis off


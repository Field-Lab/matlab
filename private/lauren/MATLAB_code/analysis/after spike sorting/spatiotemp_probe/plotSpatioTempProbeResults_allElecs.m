function plotSpatioTempProbeResults_allElecs(prePulse, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)



markerSize = 4;

offsets(offsets == 0) = [];
offsetsMs = offsets/20;

nPreElecs = length(prePulse);
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

clear temp



%%


figure('position', [100 100 800 350], 'color', [1 1 1])

aTitle = axes('position', [0.1 0.9 0.8 0.08]);
text(0.5, 0.5, ['neuron ' num2str(neuronID)], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center')
set(aTitle, 'xlim', [0 1], 'ylim', [0 1])
axis off


aLeft = axes('position', [0.08 0.15 0.42 0.7]);
aRight = axes('position', [0.54 0.15 0.42 0.7]);

%plot threshold for isolated stim pulse
axes(aLeft); hold on
plot([0 offsetsMs(end)+5.5], [threshIso threshIso], 'color', [0.8 0.8 0.8])
plot(offsetsMs(end)+5, threshIso, 'o', 'markerFaceColor', 'none', 'markerEdgeColor', [0 0 0], 'markerSize', markerSize)
plot([offsetsMs(end)+5 offsetsMs(end)+5], [threshStdIso, -1*threshStdIso] + threshIso,...
    '-', 'color', [0 0 0])
hold off

axes(aRight); hold on
plot([0 offsetsMs(end)+5.5], [threshIso threshIso], 'color', [0.8 0.8 0.8])
plot(offsetsMs(end)+5, threshIso, 'o', 'markerFaceColor', 'none', 'markerEdgeColor', [0 0 0], 'markerSize', markerSize)
plot([offsetsMs(end)+5 offsetsMs(end)+5], [threshStdIso, -1*threshStdIso] + threshIso,...
    '-', 'color', [0 0 0])
hold off

for ii = 1:nPreElecs
    nAmps = length(prePulse(ii).amps);

    for kk = 1:nAmps
        plotColorSucc = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(4,:);
        plotColorFail = (1-prePulse(ii).meanPreRespProb(kk))*plotColors(5,:);
        
        
        %plot lines connecting data points
        psThresh = prePulse(ii).threshPostSucc(:,kk);
        pfThresh = prePulse(ii).threshPostFail(:,kk);
        
        %thresholds outside measured amp range are plotted at limits of amp range
%         psThresh(isinf(psThresh) & psThresh > 0) = ampRange(2);
%         psThresh(isinf(psThresh) & psThresh < 0) = ampRange(1);
%         pfThresh(isinf(pfThresh) & pfThresh > 0) = ampRange(2);
%         pfThresh(isinf(pfThresh) & pfThresh < 0) = ampRange(1);
        
%         axes(aLeft); hold on
%         if ~any(isnan(psThresh))
%             plot(offsetsMs, psThresh, 'color', plotColorSucc)
%         else
%             for nn = 1:length(offsetsMs)-1
%                 if ~isnan(psThresh(nn)) && ~isnan(psThresh(nn)) %current and next points are valid
%                     plot(offsetsMs(nn:nn+1), psThresh(nn:nn+1), 'color', plotColorSucc)
%                 end
%             end
%         end
%         hold off
%         axes(aRight); hold on
%         if ~any(isnan(pfThresh)) && prePulse(ii).analyzedPreResponse(kk) == 1
%             plot(offsetsMs, pfThresh, 'color', plotColorFail)
%         elseif any(isnan(pfThresh))
%             for nn = 1:length(offsetsMs)-1
%                 if ~isnan(pfThresh(nn)) && ~isnan(pfThresh(nn)) %current and next points are valid
%                     plot(offsetsMs(nn:nn+1), pfThresh(nn:nn+1), 'color', plotColorFail)
%                 end
%             end
%         end
%         hold off
%         
axes(aLeft); hold on

if all(isfinite(psThresh))
    plot(offsetsMs, psThresh, 'color', plotColorSucc)
else
    for nn = 1:length(offsetsMs)-1
        if isfinite(psThresh(nn)) && isfinite(psThresh(nn)) %current and next points are valid
            plot(offsetsMs(nn:nn+1), psThresh(nn:nn+1), 'color', plotColorSucc)
        end
    end
end
hold off
axes(aRight); hold on
if all(isfinite(pfThresh)) && prePulse(ii).analyzedPreResponse(kk) == 1
    plot(offsetsMs, pfThresh, 'color', plotColorFail)
elseif ~all(isfinite(pfThresh))
    for nn = 1:length(offsetsMs)-1
        if isfinite(pfThresh(nn)) && isfinite(pfThresh(nn)) %current and next points are valid
            plot(offsetsMs(nn:nn+1), pfThresh(nn:nn+1), 'color', plotColorFail)
        end
    end
end
hold off

        
        
        for jj = 1:nOffsets
            %plots post-success thresholds
            axes(aLeft); hold on
            if prePulse(ii).threshPostSucc(jj,kk)~=0 && isfinite(prePulse(ii).threshPostSucc(jj,kk))
                plot(offsetsMs(jj), prePulse(ii).threshPostSucc(jj,kk), 'o',...
                    'markerEdgeColor', plotColorSucc, 'markerFaceColor', 'none', 'markerSize', markerSize)
                plot([offsetsMs(jj) offsetsMs(jj)],...
                    [prePulse(ii).threshStdPostSucc(jj,kk), -1*prePulse(ii).threshStdPostSucc(jj,kk)] + prePulse(ii).threshPostSucc(jj,kk),...
                    'color', plotColorSucc)
            elseif isinf(prePulse(ii).threshPostSucc(jj,kk)) %plots cases where threshold is outside measured amplitudes
                if prePulse(ii).threshPostSucc(jj,kk) > 0
                    plot(offsetsMs(jj), ampRange(2), '^', 'markerEdgeColor', plotColorSucc, 'markerSize', markerSize)
                elseif prePulse(ii).threshPostSucc(jj,kk) < 0
                    plot(offsetsMs(jj), ampRange(1), 'v', 'markerEdgeColor', plotColorSucc, 'markerSize', markerSize)
                end
            end
            hold off
            
            %plots post-failure thresholds
            axes(aRight); hold on
            if prePulse(ii).analyzedPreResponse(kk) == 1
                if prePulse(ii).threshPostFail(jj,kk)~=0 && isfinite(prePulse(ii).threshPostFail(jj,kk))
                    plot(offsetsMs(jj), prePulse(ii).threshPostFail(jj,kk), 'o',...
                        'markerEdgeColor', plotColorFail, 'markerFaceColor', 'none', 'markerSize', markerSize)
                    plot([offsetsMs(jj) offsetsMs(jj)],...
                        [prePulse(ii).threshStdPostFail(jj,kk), -1*prePulse(ii).threshStdPostFail(jj,kk)] + prePulse(ii).threshPostFail(jj,kk),...
                        'color', plotColorFail)
                elseif isinf(prePulse(ii).threshPostFail(jj,kk)) %plots cases where threshold is outside measured amplitudes
                    if prePulse(ii).threshPostFail(jj,kk) > 0
                        plot(offsetsMs(jj), ampRange(2), '^', 'markerEdgeColor', plotColorFail, 'markerSize', markerSize)
                    elseif prePulse(ii).threshPostFail(jj,kk) < 0
                        plot(offsetsMs(jj), ampRange(1), 'v', 'markerEdgeColor', plotColorFail, 'markerSize', markerSize)
                    end
                end
            end
            hold off
        end
    end   
end

axes(aLeft)
set(aLeft, 'xlim', [0 offsetsMs(end)+5.5], 'ylim', threshYLims, 'xtick', 0:2:offsetsMs(end))
ylabel('threshold (µA)')
xlabel('delay between prepulse and stim pulse (ms)')

axes(aRight)
set(aRight, 'xlim', [0 offsetsMs(end)+5.5], 'ylim', threshYLims, 'xtick', 0:2:offsetsMs(end))
ylabel('threshold (µA)')
xlabel('delay between prepulse and stim pulse (ms)')
function thresholds_full = axonStimFigurePlotter(pathToData, neuronID, patterns, pElec, sElec)

%note: this function is currently broken due to changes in analysis structure


% patterns: vector of pattern number corresponding to (in order):
%    primary electrode alone -- negative
%    primary electrode alone -- positive
%    secondary electrode alone -- negative
%    secondary electrode alone -- positive
%    primary neg. + secondary neg.
%    primary neg. + secondary pos.
%    primary pos. + secondary neg.
%    primary pos. + secondary pos.

%% for testing %%%%%%%%%%%%%

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data010';
neuronID = 559;
patterns = [49 1 2];
%patterns = [53 5 7];
pElec = 3;
sElec = 1;

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';
neuronID = 800;
%patterns = [317 249 261 297 309];
%patterns = [317 250 262 298 310];
patterns = [317 296 308];

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';
neuronID = 800;
patterns = [317 327];

pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002';
neuronID = 677;
patterns = [46 49];

neuronID = 676;
patterns = [46];

neuronID = 754;
patterns = [49];





% pElec = 35;
% sElec = 37;
% 
% preprocDataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
% pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
% neuronID = 227;
% patterns = [49:50, 57:58, 13:16];
% movies = 1:26;
% 
% load '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/old data/2008-08-27-4-data007/neuron227/latencies_complete.mat'


%% %%%%%%%%%%%%%%%%%%%%%%%%%
elecResps = cell(length(patterns), 1);
nBootstrapReps = 100;

for i = 1:length(patterns)
    temp = load([pathToData '/elecResp_n' num2str(neuronID) '_p' num2str(patterns(i))]);
    elecResp = temp.elecResp;
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', 0, 'keepLogBased', 0);
    save([pathToData 'elecResp_n' num2str(neuronID) '_p' num2str(patterns(i))], 'elecResp')
    elecResps{i} = elecResp;
end

[xCoords yCoords] = getElectrodeCoords61();
ei = elecResps{1}.cells.mainEI;

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for i = 1:64
    if ~(i==9||i==25||i==57)
        eiAmps(i) = max(max(abs(ei(i,:))));
    end
end

eiAmpsNormalized = eiAmps/max(eiAmps);

%% plotting
% figure
% set(gcf, 'color', 'white')

% primColor = [0 0 1];
% secColor = [1 0 0];
% sameColor = [1 0 1];
% oppColor = [0 0.8 0.5];
% 
% primColor = [74 138 247]/255;

% vintage candy shop
% primColor = [100/255 168/255 189/255];
% sameColor = [135/255 184/255 67/255];
% oppColor = [227/255 100/255 54/255];

%primColor = [101 52 255]/255;% purple
%sameColor = [52 198 247]/255; %aqua
%oppColor = [247 52 198]/255; %calm magenta

primColor = [50 70 247]/255; %blue
secColor = [.8 .05 0.05];  %rust
sameColor = [90 156 0]/255; %pale grass
oppColor = [255 124 59]/255;  %salmon

plotColors = [primColor; sameColor; oppColor];

% plotColors(1, :) = [0 0 0];
% plotColors(4, :) = [1 0.4 0.4]; %pink
% plotColors(2, :) = [0.8 0 0]; %dark red
% plotColors(3, :) = [0 0 0.8]; %dark blue
% plotColors(5, :) = [0.4 0.4 1]; %light blue

%% ei portion

% subplot(2,1,1)
% hold on
% 
% for i = 1:64
%     if ~(i==9||i==25||i==57)
%         if eiAmpsNormalized(i)>.02
%             plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
%         end
%     end
% end
% 
% %plot(xCoords(pElec), yCoords(pElec), 'o', 'MarkerFaceColor', primColor, 'MarkerEdgeColor',
% %primColor)
% %plot(xCoords(sElec), yCoords(sElec), 'o', 'MarkerFaceColor', secColor, 'MarkerEdgeColor', secColor)
% %plot(xCoords(sElec)+0.4, yCoords(sElec), 'o', 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', [1 0 1])
% 
% axis equal
% set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
% axis off
% hold off



%% stimulation legend
% subplot(2,2,2)
% hold on
% 
% markerPositions = linspace(8, 2, 6);
% 
% % plotColors = [0 0 1;
% %     0 1 0.5;
% %     1 0 1;
% %     0 1 1;
% %     1 0.5 0];
% 
% % for i = 1:5
% %     plot(1, markerPositions(i), 's', 'MarkerFaceColor', plotColors(i,:), 'MarkerEdgeColor', plotColors(i,:))
% % end
% 
% text(1, 10, ['neuron ', num2str(neuronID), ', electrode pair ', num2str(pElec), ', ', num2str(sElec)])
% 
% plot(1, markerPositions(1), 's', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
% text(2, markerPositions(1), 'primary electrode')
% 
% plot(1, markerPositions(2), 's', 'MarkerFaceColor', secColor, 'MarkerEdgeColor', secColor)
% text(2, markerPositions(2), 'secondary electrode')
% 
% plot(1, markerPositions(4), 's', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
% text(2, markerPositions(4), 'primary electrode alone')
% 
% plot(1, markerPositions(5), 's', 'MarkerFaceColor', sameColor, 'MarkerEdgeColor', sameColor)
% text(2, markerPositions(5), 'primary and secondary, same polarity')
% 
% plot(1, markerPositions(6), 's', 'MarkerFaceColor', oppColor, 'MarkerEdgeColor', oppColor)
% text(2, markerPositions(6), 'primary and secondary, opposite polarity')
%     
% 
% set(gca, 'xlim', [0 10], 'ylim', [0 10])
% axis off
% hold off

%% response curves

thresholds_full = zeros(8,1);

figure('position', [100 100 161.8*2 200])
axes('position', [0.2 0.2 0.7 0.7])
hold on


% 1   primary electrode alone -- negative
% 5   primary neg. + secondary neg.
% 6   primary neg. + secondary pos.

%plotColors = [primColor; sameColor; oppColor];

erfFitPlotter(elecResps, gca, plotColors, 1)

% plots tick marks at threshold locations
hold on
for i = 1:length(elecResps)
    plot([elecResps{i}.analysis.threshold elecResps{i}.analysis.threshold], [0 0.1], 'Color', plotColors(i,:))
end
hold off

%set(gca, 'Xlim', [0 1], 'Ylim', [0 1], 'XTick', [0 0.2 0.4 0.6 0.8 1], 'YTick', [0 0.5 1], 'YTickLabel', [0 50 100])
%set(gca, 'Xlim', [0.45 1.8], 'Ylim', [0 1], 'YTick', [0 0.5 1])
set(gca, 'position', [0.2 0.2 0.7 0.7], 'Xlim', [0.1 1.2], 'Ylim', [0 1], 'YTick', [0 0.5 1], 'XScale', 'linear')
xlabel('current amplitude (\muA)')
ylabel('spike probability')

% text(thresholds(1), 1.1, '-', 'fontsize', 16)
% text(thresholds(2), 1.1, '-', 'fontsize', 16)
% text(thresholds(3), 1.1, '+', 'fontsize', 16)

% thresholds_full(1) = thresholds(1);
% thresholds_full(5) = thresholds(2);
% thresholds_full(6) = thresholds(3);


hold off

% subplot(2,2,4)
% hold on

% posPrimPatterns = patterns([2 8 7]);
% if isLatencies
%     thresholds = erfFitPlotter(posPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor], erfFitDetails);
% else
%     thresholds = erfFitPlotter(posPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor], erfFitDetails);
% end
% axis tight
% %set(gca, 'Xlim', [0 1], 'Ylim', [0 1], 'XTick', [0 0.2 0.4 0.6 0.8 1], 'YTick', [0 0.5 1], 'YTickLabel', [0 50 100])
% set(gca, 'Ylim', [0 1], 'YTick', [0 0.5 1], 'YTickLabel', [0 0.5 1], 'fontsize', 18)



% 
% text(thresholds(1), 1.1, '+', 'fontsize', 16)
% text(thresholds(2), 1.1, '+', 'fontsize', 16)
% text(thresholds(3), 1.1, '-', 'fontsize', 16)
% 
% thresholds_full(2) = thresholds(1);
% thresholds_full(8) = thresholds(2);
% thresholds_full(7) = thresholds(3);

% hold off


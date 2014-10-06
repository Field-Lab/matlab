function thresholds_full = somaStimSpecialFigurePlotter(responses, pElec, sElec, patterns, movies, preprocDataPath, pathToEi, neuronID, erfFitDetails)

% patterns: vector of pattern number corresponding to (in order):
%    primary electrode alone -- negative
%    primary electrode alone -- positive
%    secondary electrode alone -- negative
%    secondary electrode alone -- positive
%    primary neg. + secondary neg.
%    primary neg. + secondary pos.
%    primary pos. + secondary neg.
%    primary pos. + secondary pos.


% a temporary cloog
if exist('erfFitDetails', 'var')
    isLatencies = 0;
else
    isLatencies = 1;
    %latencies = responses;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%

[xCoords yCoords] = getElectrodeCoords61();

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
clear eiFile

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for i = 1:64
    if ~(i==9||i==25||i==57)
        eiAmps(i) = max(max(abs(ei(1,i+1,:))));
    end
end

eiAmpsNormalized = eiAmps/max(eiAmps);

%% plotting
figure
set(gcf, 'color', 'white')


primColor = [50 70 247]/255; %blue
secColor = [.8 .05 0.05];  %rust
sameColor = [90 156 0]/255; %pale grass
oppColor = [255 124 59]/255;  %salmon

%% ei portion
subplot(2,2,1)
hold on

for i = 1:64
    if ~(i==9||i==25||i==57)
        if eiAmpsNormalized(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
        end
    end
end

plot(xCoords(pElec), yCoords(pElec), 'o', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
plot(xCoords(sElec), yCoords(sElec), 'o', 'MarkerFaceColor', secColor, 'MarkerEdgeColor', secColor)
%plot(xCoords(sElec)+0.4, yCoords(sElec), 'o', 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', [1 0 1])


set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off



%% stimulation legend
subplot(2,2,2)
hold on

markerPositions = linspace(8, 2, 6);

% plotColors = [0 0 1;
%     0 1 0.5;
%     1 0 1;
%     0 1 1;
%     1 0.5 0];

% for i = 1:5
%     plot(1, markerPositions(i), 's', 'MarkerFaceColor', plotColors(i,:), 'MarkerEdgeColor', plotColors(i,:))
% end

text(1, 10, ['neuron ', num2str(neuronID), ', electrode pair ', num2str(pElec), ', ', num2str(sElec)])

plot(1, markerPositions(1), 's', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
text(2, markerPositions(1), 'primary electrode')

plot(1, markerPositions(2), 's', 'MarkerFaceColor', secColor, 'MarkerEdgeColor', secColor)
text(2, markerPositions(2), 'secondary electrode')

plot(1, markerPositions(4), 's', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
text(2, markerPositions(4), 'primary electrode alone')

plot(1, markerPositions(5), 's', 'MarkerFaceColor', sameColor, 'MarkerEdgeColor', sameColor)
text(2, markerPositions(5), 'primary and secondary, same polarity')

plot(1, markerPositions(6), 's', 'MarkerFaceColor', oppColor, 'MarkerEdgeColor', oppColor)
text(2, markerPositions(6), 'primary and secondary, opposite polarity')
    

set(gca, 'xlim', [0 10], 'ylim', [0 10])
axis off
hold off

%% response curves

thresholds_full = zeros(8,1);

subplot(2,2,3)
hold on

% 1   primary electrode alone -- negative
% 2   primary electrode alone -- positive
% 3   secondary electrode alone -- negative
% 4   secondary electrode alone -- positive
% 5   primary neg. + secondary neg.
% 6   primary neg. + secondary pos.
% 7   primary pos. + secondary neg.
% 8   primary pos. + secondary pos.

negPrimPatterns = patterns([1 5 6]);
if isLatencies
    thresholds = erfFitPlotter(negPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor]);
else
    thresholds = erfFitPlotter(negPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor], erfFitDetails);
end
axis tight

text(thresholds(1), 1.1, '-')
text(thresholds(2), 1.1, '-')
text(thresholds(3), 1.1, '+')

thresholds_full(1) = thresholds(1);
thresholds_full(5) = thresholds(2);
thresholds_full(6) = thresholds(3);

hold off

subplot(2,2,4)
hold on

posPrimPatterns = patterns([2 8 7]);
if isLatencies
    thresholds = erfFitPlotter(posPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor], erfFitDetails);
else
    thresholds = erfFitPlotter(posPrimPatterns, movies, preprocDataPath, responses, [primColor; sameColor; oppColor], erfFitDetails);
end
axis tight

text(thresholds(1), 1.1, '+')
text(thresholds(2), 1.1, '+')
text(thresholds(3), 1.1, '-')

thresholds_full(2) = thresholds(1);
thresholds_full(8) = thresholds(2);
thresholds_full(7) = thresholds(3);

hold off


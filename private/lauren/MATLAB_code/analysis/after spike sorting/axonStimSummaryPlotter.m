%function axonStimSummaryPlotter(latencies, patterns, movies, preprocDataPath, pathToEi, neuronID)

%% for testing %%%%%%%%%%%%%

preprocDataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 227;
patterns = [1:62];
movies = 1:26;


%% %%%%%%%%%%%%%%%%%%%%%%%%%


[xCoords yCoords] = getElectrodeCoords61();

% get the electrodes that are active in each pattern and the polarity of their stim
elecInPattern = cell(length(patterns), 1);
polarities = cell(length(patterns), 1);
for i = 1:length(patterns)
    [amps elecInPattern{i}] = getStimAmps(preprocDataPath, patterns(i), 1);
    polarities{i} = amps./abs(amps);
end

% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
clear eiFile

% calculate maximum waveform value on each electrode (absolute value)
eiAmps = zeros(64, 1);
for i = 1:64
    if ~(i==9||i==25||i==57)
        eiAmps(i) = max(max(ei(1,i+1,:)));
    end
end

eiAmpsNormalized = eiAmps/max(eiAmps);

%% plotting
figure
set(gcf, 'color', 'white')
plotColors = jet(length(patterns));


%% ei portion
subplot(2,2,1)
hold on

for i = 1:64
    if ~(i==9||i==25||i==57)
        plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
    end
end

for i = 1:length(patterns)
    elec = elecInPattern{i};
    for j = 1:length(elec)
        text(xCoords(elec(j))-0.3, yCoords(elec(j))-0.1, num2str(j), 'color', 'w')
    end
end


axis equal
axis off
hold off

%% stimulation legend
subplot(2,2,2)
hold on

markerPositions = linspace(2, 8, length(patterns));


for i = 1:length(patterns)
    elec = elecInPattern{i};
    polar = polarities{i};
    
    plot(1, markerPositions(i), 's', 'MarkerFaceColor', plotColors(i,:), 'MarkerEdgeColor', plotColors(i,:))
    for j = 1:length(polar)
        text(2*j, markerPositions(i), num2str(j))
        if polar(j)==1
            text(1+2*j, markerPositions(i), '+')
        else
            text(1+2*j, markerPositions(i), '-')
        end
    end
    text(4, markerPositions(i), '')
end

set(gca, 'xlim', [0 10], 'ylim', [0 10])
axis off
hold off

%% response curves
subplot(2,1,2)
hold on

erfFitPlotter(patterns, movies, preprocDataPath, latencies, plotColors)
axis tight

hold off


function maxEiAmp = somaStimFigurePlotter(responses, patternNos, pathToEi, neuronID, stimAmps, projection)


%% for testing %%%%%%%%%%%%%


% preprocDataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
% pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
% neuronID = 227;
% patterns = [49:50, 57:58, 13:16];
% movies = 1:26;
% 
% load '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/old data/2008-08-27-4-data007/neuron227/latencies_complete.mat'


%% %%%%%%%%%%%%%%%%%%%%%%%%%

nPatterns = length(patternNos);
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

maxEiAmp = max(eiAmps);
eiAmpsNormalized = eiAmps/maxEiAmp;


%% plotting
figure('Position',[100 100 300 400])
set(gcf, 'color', 'white')

% vintage candy shop
% primColor = [100/255 168/255 189/255];
% sameColor = [135/255 184/255 67/255];
% oppColor = [227/255 100/255 54/255];

color{1} = [50 70 247]/255; %blue
color{2} = [.8 .05 0.05];  %rust
color{3} = [90 156 0]/255; %pale grass
color{4} = [255 124 59]/255; %salmon
color{5} = [101 52 255]/255; %purple
color{6} = [52 198 247]/255; %aqua
color{7} = [52 198 247]/255; %calm magenta

%% ei portion
subplot(2,1,1)
hold on

for i = 1:64
    if ~(i==9||i==25||i==57)
        if eiAmpsNormalized(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
        end
    end
end

for i = 1:nPatterns
    pattern = patternNos(i);
    plot(xCoords(pattern), yCoords(pattern), 'o', 'MarkerFaceColor', color{i}, 'MarkerEdgeColor', color{i})
end

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis equal
axis off

title(['neuron ', num2str(neuronID)])

hold off



%% response curves

subplot(2,1,2)
hold on

for i = 1:nPatterns
    plot(stimAmps, responses,'o','MarkerEdgeColor', color{i},'MarkerFaceColor', color{i})
    current = plot(projection(1,:), projection(2,:));
    set(findobj(current,'Type','line'),'Color', color{i})
end
hold off
xlabel('current amplitude (\muA)')
ylabel('response rate')

axis tight

hold off




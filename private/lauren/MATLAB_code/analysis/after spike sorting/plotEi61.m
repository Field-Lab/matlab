function [eiAmps scaleFactor] = plotEi61(pathToEi, neuronIDs, varargin)

% plots a neurons EI on the 61-array layout in specified axes, marking any electrodes in markElecs
% with a colored dot

% if more than one neuron is specified, plots each neuron's ei as a colored outline rather than
% solid black dots
%

p = inputParser;

color(1,:) = [50 70 247]/255; %blue
color(2,:) = [.8 .05 0.05];  %rust
color(3,:) = [90 156 0]/255; %pale grass
color(4,:) = [255 124 59]/255; %salmon
color(5,:) = [101 52 255]/255; %purple
color(6,:) = [52 198 247]/255; %aqua
color(7,:) = [238 55 128]/255; %calm magenta



p.addRequired('pathToEi', @ischar)
p.addRequired('neuronIDs', @isnumeric)

p.addParamValue('markElecs', [], @isnumeric)
p.addParamValue('axesH', gca, @ishandle)
p.addParamValue('neuronColors', color, @isnumeric) %only works if more than 1 neuron specified
p.addParamValue('markElecColors', color, @isnumeric)
p.addParamValue('universalScale', 0, @(x)any(x==[0 1])) %for multiple neurons, scale all ei amplitudes to universal max
p.addParamValue('eiThresh', .02, @isnumeric) %minimum ei amplitude displayed on a given electrode (relative to max ei value)
p.addParamValue('manualScale', [], @isnumeric)
p.addParamValue('eiData', {}, @iscell) %if supplied, this data is used instead of data from .ei file (must be same format as eiAmps)

p.parse(pathToEi, neuronIDs, varargin{:})

markElecs = p.Results.markElecs;
axesH = p.Results.axesH;
neuronColors = p.Results.neuronColors;
markElecColors = p.Results.markElecColors;
universalScale = p.Results.universalScale;
eiThresh = p.Results.eiThresh;
manualScale = p.Results.manualScale;
eiData = p.Results.eiData;






[xCoords yCoords] = getElectrodeCoords61();

if isempty(eiData)
    nNeuron = length(neuronIDs);
    % get neuron's average spike waveforms from VISION analysis
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
    ei = cell(nNeuron, 1);
    for i = 1:nNeuron
        ei{i} = eiFile.getImage(neuronIDs(i)); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
    end
    clear eiFile

    % calculate maximum waveform value on each electrode (absolute value)

    eiAmps = cell(nNeuron, 1);
    maxEiAmps = zeros(nNeuron, 1);
    for j = 1:nNeuron
        eiAmps{j} = zeros(64, 1);
        for i = 1:64
            if ~(i==9||i==25||i==57)
                eiAmps{j}(i) = max(max(abs(ei{j}(1,i+1,:))));
            end
        end
        maxEiAmps(j) = max(eiAmps{j});
    end
else
    nNeuron = length(eiData);
    maxEiAmps = zeros(nNeuron, 1);
    eiAmps = eiData;
    for j = 1:nNeuron
        maxEiAmps(j) = max(eiAmps{j});
    end
end

if isempty(manualScale)
    maxEiAmp = max(maxEiAmps);
    scaleFactor = 1/maxEiAmp;
else
    scaleFactor = manualScale;
end

eiAmpsNormalized = cell(nNeuron, 1);
for j = 1:nNeuron
    if universalScale || ~isempty(manualScale)
        eiAmpsNormalized{j} = eiAmps{j}*scaleFactor;
    else
        eiAmpsNormalized{j} = eiAmps{j}/maxEiAmps(j);
    end
end

axes(axesH)
hold on

for j = 1:nNeuron
    for i = 1:64
        if ~(i==9||i==25||i==57)
            if eiAmpsNormalized{j}(i)>eiThresh
                if nNeuron > 1 
                    plot(xCoords(i), yCoords(i), 'o','MarkerSize', ceil(eiAmpsNormalized{j}(i)*20), 'MarkerEdgeColor', neuronColors(j,:))
                else
                    plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized{j}(i)*20), 'MarkerFaceColor', 'k')
                end
            end
        end
    end
end

if exist('markElecs','var')
    for i = 1:length(markElecs)
        plot(xCoords(markElecs(i)), yCoords(markElecs(i)), 'o', 'MarkerFaceColor', color(i,:), 'MarkerEdgeColor', markElecColors(i,:))
    end
end

% array outline
plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')

set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
axis equal
axis off

hold off
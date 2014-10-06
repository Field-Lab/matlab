blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon

errordlg('This script needs to be updated to work with the new erfFitter (now uses maximum likelihood fitting)')

n31responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_31.csv');
n557responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_557.csv');
n801responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_801.csv');
n77responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_77.csv');
n316responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_316.csv');
n801_2responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_801-2.csv');
n558responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/somal_stim/soma_responses_558.csv');


nstimAmpsFull = [0.5055 0.5561 0.6066 0.6572 0.7330 0.8088 0.8846 0.9857 1.1043 1.2047 1.3051 1.4055...
    1.6063 1.7067 1.9075 2.1083 2.3091 2.5098 2.8110 3.0118 3.4134 3.7146 4.1161]; % 11-10-3
n801stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043]; % 8-27-4
n557stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043]; % 8-27-4
n31stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043...
    1.2047 1.3051 1.4055 1.6063 1.7067 1.9075 2.1083]; % 8-26-0

n77stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857];%8-27-4 cluster
n316stimAmpsFull = n77stimAmpsFull;
n801_2stimAmpsFull = n77stimAmpsFull;
n558stimAmpsFull = n77stimAmpsFull;

n31patternNos = n31responsesFull(:,1);
n31responses = cell(length(n31patternNos), 1);
n31stimAmps = cell(length(n31patternNos), 1);
n31params = cell(length(n31patternNos), 1);
n31projections = cell(length(n31patternNos), 1);
n31thresh = zeros(length(n31patternNos), 1);

n557patternNos = n557responsesFull(:,1);
n557responses = cell(length(n557patternNos), 1);
n557stimAmps = cell(length(n557patternNos), 1);
n557params = cell(length(n557patternNos), 1);
n557projections = cell(length(n557patternNos), 1);
n557thresh = zeros(length(n557patternNos), 1);

n801patternNos = n801responsesFull(:,1);
n801responses = cell(length(n801patternNos), 1);
n801stimAmps = cell(length(n801patternNos), 1);
n801params = cell(length(n801patternNos), 1);
n801projections = cell(length(n801patternNos), 1);
n801thresh = zeros(length(n801patternNos), 1);

n77patternNos = n77responsesFull(:,1);
n77responses = cell(length(n77patternNos), 1);
n77stimAmps = cell(length(n77patternNos), 1);
n77params = cell(length(n77patternNos), 1);
n77projections = cell(length(n77patternNos), 1);
n77thresh = zeros(length(n77patternNos), 1);

n316patternNos = n316responsesFull(:,1);
n316responses = cell(length(n316patternNos), 1);
n316stimAmps = cell(length(n316patternNos), 1);
n316params = cell(length(n316patternNos), 1);
n316projections = cell(length(n316patternNos), 1);
n316thresh = zeros(length(n316patternNos), 1);

n801_2patternNos = n801_2responsesFull(:,1);
n801_2responses = cell(length(n801_2patternNos), 1);
n801_2stimAmps = cell(length(n801_2patternNos), 1);
n801_2params = cell(length(n801_2patternNos), 1);
n801_2projections = cell(length(n801_2patternNos), 1);
n801_2thresh = zeros(length(n801_2patternNos), 1);

n558patternNos = n558responsesFull(:,1);
n558responses = cell(length(n558patternNos), 1);
n558stimAmps = cell(length(n558patternNos), 1);
n558params = cell(length(n558patternNos), 1);
n558projections = cell(length(n558patternNos), 1);
n558thresh = zeros(length(n558patternNos), 1);

for i = 1:length(n31patternNos)
    n31responses{i} = n31responsesFull(i, 2:end);
    n31stimAmps{i} = n31stimAmpsFull;
    for j = length(n31responses{i}):-1:1
        if n31responses{i}(j) > 1.1
            n31responses{i}(j) = [];
            n31stimAmps{i}(j) = [];
        end
    end
    [n31params{i} n31projections{i}] = erfFitter([n31stimAmps{i}; n31responses{i}], 1, 1);
    n31thresh(i) = -n31params{i}(2)/n31params{i}(1);
end

for i = 1:length(n557patternNos)
    n557responses{i} = n557responsesFull(i, 2:end);
    n557stimAmps{i} = n557stimAmpsFull;
    for j = length(n557responses{i}):-1:1
        if n557responses{i}(j) > 1.1
            n557responses{i}(j) = [];
            n557stimAmps{i}(j) = [];
        end
    end
    [n557params{i} n557projections{i}] = erfFitter([n557stimAmps{i}; n557responses{i}], 1, 1);
    n557thresh(i) = -n557params{i}(2)/n557params{i}(1);
end

for i = 1:length(n801patternNos)
    n801responses{i} = n801responsesFull(i, 2:end);
    n801stimAmps{i} = n801stimAmpsFull;
    for j = length(n801responses{i}):-1:1
        if n801responses{i}(j) > 1.1
            n801responses{i}(j) = [];
            n801stimAmps{i}(j) = [];
        end
    end
    [n801params{i} n801projections{i}] = erfFitter([n801stimAmps{i}; n801responses{i}], 1, 1);
    n801thresh(i) = -n801params{i}(2)/n801params{i}(1);
end

for i = 1:length(n77patternNos)
    n77responses{i} = n77responsesFull(i, 2:end);
    n77stimAmps{i} = n77stimAmpsFull;
    for j = length(n77responses{i}):-1:1
        if n77responses{i}(j) > 1.1
            n77responses{i}(j) = [];
            n77stimAmps{i}(j) = [];
        end
    end
    [n77params{i} n77projections{i}] = erfFitter([n77stimAmps{i}; n77responses{i}], 1, 1);
    n77thresh(i) = -n77params{i}(2)/n77params{i}(1);
end

for i = 1:length(n316patternNos)
    n316responses{i} = n316responsesFull(i, 2:end);
    n316stimAmps{i} = n316stimAmpsFull;
    for j = length(n316responses{i}):-1:1
        if n316responses{i}(j) > 1.1
            n316responses{i}(j) = [];
            n316stimAmps{i}(j) = [];
        end
    end
    [n316params{i} n316projections{i}] = erfFitter([n316stimAmps{i}; n316responses{i}], 1, 1);
    n316thresh(i) = -n316params{i}(2)/n316params{i}(1);
end

for i = 1:length(n801_2patternNos)
    n801_2responses{i} = n801_2responsesFull(i, 2:end);
    n801_2stimAmps{i} = n801_2stimAmpsFull;
    for j = length(n801_2responses{i}):-1:1
        if n801_2responses{i}(j) > 1.1
            n801_2responses{i}(j) = [];
            n801_2stimAmps{i}(j) = [];
        end
    end
    [n801_2params{i} n801_2projections{i}] = erfFitter([n801_2stimAmps{i}; n801_2responses{i}], 1, 1);
    n801_2thresh(i) = -n801_2params{i}(2)/n801_2params{i}(1);
end


for i = 1:length(n558patternNos)
    n558responses{i} = n558responsesFull(i, 2:end);
    n558stimAmps{i} = n558stimAmpsFull;
    for j = length(n558responses{i}):-1:1
        if n558responses{i}(j) > 1.1
            n558responses{i}(j) = [];
            n558stimAmps{i}(j) = [];
        end
    end
    [n558params{i} n558projections{i}] = erfFitter([n558stimAmps{i}; n558responses{i}], 1, 1);
    n558thresh(i) = -n558params{i}(2)/n558params{i}(1);
end

% patterns: vector of pattern number corresponding to (in order):
%    primary electrode alone -- negative
%    primary electrode alone -- positive
%    secondary electrode alone -- negative
%    secondary electrode alone -- positive
%    primary neg. + secondary neg.
%    primary neg. + secondary pos.
%    primary pos. + secondary neg.
%    primary pos. + secondary pos.

%% plotting: neuron 31
movies = []; %for now, a required argument to plot, but not used

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data010';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data009/data009.ei';
neuronID = 31;
responses = n31responses;
erfFitDetails.stimAmps    = n31stimAmps;
erfFitDetails.projections = n31projections;
erfFitDetails.thresholds  = n31thresh;
erfFitDetails.patternNos  = n31patternNos;

nPairs = 1;
n31threshFull = cell(nPairs, 1);
n31pairIDs = cell(nPairs, 1);

% electrode pair 3, 5 (primary = 3)
pElec = 3;
sElec = 5;
patterns = [49 50 55 56 9:12];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n31threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

%% plotting: neuron 557

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 557;
responses = n557responses;
erfFitDetails.stimAmps    = n557stimAmps;
erfFitDetails.projections = n557projections;
erfFitDetails.thresholds  = n557thresh;
erfFitDetails.patternNos  = n557patternNos;

nPairs = 1;
n557threshFull = cell(nPairs, 1);
n557pairIDs = cell(nPairs, 1);

% electrode pair 38, 39 (primary = 38)
pElec = 38;
sElec = 39;
patterns = [59:62 45:48];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n557threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

%% plotting: neuron 801

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 801;
responses = n801responses;
erfFitDetails.stimAmps    = n801stimAmps;
erfFitDetails.projections = n801projections;
erfFitDetails.thresholds  = n801thresh;
erfFitDetails.patternNos  = n801patternNos;

nPairs = 1;
n801threshFull = cell(nPairs, 1);
n801pairIDs = cell(nPairs, 1);

% electrode pair 49, 46 (primary = 49)
pElec = 49;
sElec = 46;
patterns = [121:122 111:112 79 81 80 82];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n801threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

%% plotting: neuron 77

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data005';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 77;
responses = n77responses;
erfFitDetails.stimAmps    = n77stimAmps;
erfFitDetails.projections = n77projections;
erfFitDetails.thresholds  = n77thresh;
erfFitDetails.patternNos  = n77patternNos;

nPairs = 1;
n77threshFull = cell(nPairs, 1);
n77pairIDs = cell(nPairs, 1);

% electrode pair 1, 3 (primary = 1)
pElec = 1;
sElec = 3;
patterns = [99 100 97 98 73 79 85 91];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n77threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end



%% plotting: neuron 316

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data005';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 316;
responses = n316responses;
erfFitDetails.stimAmps    = n316stimAmps;
erfFitDetails.projections = n316projections;
erfFitDetails.thresholds  = n316thresh;
erfFitDetails.patternNos  = n316patternNos;

nPairs = 2;
n316threshFull = cell(nPairs, 1);
n316pairIDs = cell(nPairs, 1);

% electrode pair 22, 24 (primary = 24)
pElec = 24;
sElec = 22;
patterns = [215 216 207 208 186 192 198 204];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n316threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

% electrode pair 22, 28 (primary = 28)
pElec = 28;
sElec = 22;
patterns = [219 220 207 208 188 194 200 206];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n316threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
end

%% plotting: neuron 801-2 (from cluster stim)

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data005';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 801;
responses = n801_2responses;
erfFitDetails.stimAmps    = n801_2stimAmps;
erfFitDetails.projections = n801_2projections;
erfFitDetails.thresholds  = n801_2thresh;
erfFitDetails.patternNos  = n801_2patternNos;

nPairs = 2;
n801_2threshFull = cell(nPairs, 1);
n801_2pairIDs = cell(nPairs, 1);

% electrode pair 49, 46 (primary = 49)
pElec = 49;
sElec = 46;
patterns = [317 318 321 322 294 306 300 315];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n801_2threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

% electrode pair 49, 53 (primary = 49)
pElec = 49;
sElec = 53;
patterns = [317 318 327 328 297 309 303 315];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n801_2threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
end
n801_2threshFull{2}(7) = 1.1043;



%% plotting: neuron 558 (from cluster stim)

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data008';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data007/data007.ei';
neuronID = 558;
responses = n558responses;
erfFitDetails.stimAmps    = n558stimAmps;
erfFitDetails.projections = n558projections;
erfFitDetails.thresholds  = n558thresh;
erfFitDetails.patternNos  = n558patternNos;

nPairs = 6;
n558threshFull = cell(nPairs, 1);
n558pairIDs = cell(nPairs, 1);

% electrode pair 44, 36 (primary = 44)
pElec = 44;
sElec = 36;
patterns = [317 318 319 320 293 305 299 311];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end

% electrode pair 44, 38 (primary = 44)
pElec = 44;
sElec = 38;
patterns = [317 318 321 322 294 306 300 312];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
end

% electrode pair 44, 43 (primary = 44)
pElec = 44;
sElec = 43;
patterns = [317 318 323 324 295 307 301 313];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{3}(i) = erfFitDetails.thresholds(patternIndex);
end

% electrode pair 44, 46 (primary = 44)
pElec = 44;
sElec = 46;
patterns = [317 318 325 326 296 308 302 314];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{4}(i) = erfFitDetails.thresholds(patternIndex);
end


% electrode pair 44, 47 (primary = 44)
pElec = 44;
sElec = 47;
patterns = [317 318 327 328 297 309 303 315];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{5}(i) = erfFitDetails.thresholds(patternIndex);
end


% electrode pair 44, 49 (primary = 44)
pElec = 44;
sElec = 49;
patterns = [317 318 329 330 298 310 304 316];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{6}(i) = erfFitDetails.thresholds(patternIndex);
end

%% special plot for neuron 558

pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data007/data007.ei';
neuronID = 558;
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

figure
subplot(2,2,1) %to get correct proportions
hold on

for i = 1:64
    if ~(i==9||i==25||i==57)
        if eiAmpsNormalized(i)>.02
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
        end
    end
end

plot(xCoords(44), yCoords(44), 'o', 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
plot(xCoords(36), yCoords(36), 'o', 'MarkerEdgeColor', rust)
plot(xCoords(38), yCoords(38), 'o', 'MarkerEdgeColor', rust)
plot(xCoords(43), yCoords(43), 'o', 'MarkerEdgeColor', rust)
plot(xCoords(46), yCoords(46), 'o', 'MarkerEdgeColor', rust)
plot(xCoords(47), yCoords(47), 'o', 'MarkerEdgeColor', rust)
plot(xCoords(49), yCoords(49), 'o', 'MarkerEdgeColor', rust)

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8], 'box', 'on')
axis equal
axis off
%set(gca, 'Xtick', [], 'Ytick', [])
hold off





%% defining electrode angles (vertex = primary, angle drawn from positive x axis as plotted)
n77pairAngles = 30;
n801pairAngles = 210;
n557pairAngles = 210;
n316pairAngles = [90 30];
n31pairAngles = 30;
n801_2pairAngles = [210 90];
n558pairAngles = [330 270 210 150 30 90];

%% defining angle of axon (vertex = soma, angle draw from positive x axis as plotted)
n77axonAngle = 187.7; %average of n227, n243 from same piece
n801axonAngle = 187.7;
n557axonAngle = 187.7;
n316axonAngle = 187.7;
n31axonAngle = 31.3; % n558, from same piece
n801_2axonAngle = 187.7;
n558axonAngle = 31.3;
%n227axonAngle = 186.6;
%n243axonAngle = 188.8;
% n334axonAngle = 109.3;
% n558axonAngle = 31.3;
% n886axonAngle = 111.2;

n77diffAngles = abs(n77pairAngles - ones(size(n77pairAngles))*n77axonAngle);
n801diffAngles = abs(n801pairAngles - ones(size(n801pairAngles))*n801axonAngle);
n557diffAngles = abs(n557pairAngles - ones(size(n557pairAngles))*n557axonAngle);
n316diffAngles = abs(n316pairAngles - ones(size(n316pairAngles))*n316axonAngle);
n31diffAngles = abs(n31pairAngles - ones(size(n31pairAngles))*n31axonAngle);
n801_2diffAngles = abs(n801_2pairAngles - ones(size(n801_2pairAngles))*n801_2axonAngle);
n558diffAngles = abs(n558pairAngles - ones(size(n558pairAngles))*n558axonAngle);

%% converts angles to equivalent within (1, 90)

for i = 1:length(n77diffAngles)
    if n77diffAngles(i) > 180;
        n77diffAngles(i) = 360 - n77diffAngles(i);
    end
    if n77diffAngles(i) > 90
        n77diffAngles(i) = 180 - n77diffAngles(i);
    end
end

for i = 1:length(n801diffAngles)
    if n801diffAngles(i) > 180;
        n801diffAngles(i) = 360 - n801diffAngles(i);
    end
    if n801diffAngles(i) > 90
        n801diffAngles(i) = 180 - n801diffAngles(i);
    end
end

for i = 1:length(n557diffAngles)
    if n557diffAngles(i) > 180;
        n557diffAngles(i) = 360 - n557diffAngles(i);
    end
    if n557diffAngles(i) > 90
        n557diffAngles(i) = 180 - n557diffAngles(i);
    end
end

for i = 1:length(n316diffAngles)
    if n316diffAngles(i) > 180;
        n316diffAngles(i) = 360 - n316diffAngles(i);
    end
    if n316diffAngles(i) > 90
        n316diffAngles(i) = 180 - n316diffAngles(i);
    end
end

for i = 1:length(n31diffAngles)
    if n31diffAngles(i) > 180;
        n31diffAngles(i) = 360 - n31diffAngles(i);
    end
    if n31diffAngles(i) > 90
        n31diffAngles(i) = 180 - n31diffAngles(i);
    end
end

for i = 1:length(n801_2diffAngles)
    if n801_2diffAngles(i) > 180;
        n801_2diffAngles(i) = 360 - n801_2diffAngles(i);
    end
    if n801_2diffAngles(i) > 90
        n801_2diffAngles(i) = 180 - n801_2diffAngles(i);
    end
end

for i = 1:length(n558diffAngles)
    if n558diffAngles(i) > 180;
        n558diffAngles(i) = 360 - n558diffAngles(i);
    end
    if n558diffAngles(i) > 90
        n558diffAngles(i) = 180 - n558diffAngles(i);
    end
end

%% summary plot

PNegSNegPoints = [];
PNegSPosPoints = [];
PPosSNegPoints = [];
PPosSPosPoints = [];

for i = 1:length(n77diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n77diffAngles(i) n77threshFull{i}(5)/n77threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n77diffAngles(i) n77threshFull{i}(6)/n77threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n77diffAngles(i) n77threshFull{i}(8)/n77threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n77diffAngles(i) n77threshFull{i}(7)/n77threshFull{i}(2)]; %#ok<AGROW>
end

% for i = 1:length(n801diffAngles)
%    PNegSNegPoints = [PNegSNegPoints; n801diffAngles(i) n801threshFull{i}(5)/n801threshFull{i}(1)]; %#ok<AGROW>
%    PNegSPosPoints = [PNegSPosPoints; n801diffAngles(i) n801threshFull{i}(6)/n801threshFull{i}(1)]; %#ok<AGROW>
%    PPosSPosPoints = [PPosSPosPoints; n801diffAngles(i) n801threshFull{i}(8)/n801threshFull{i}(2)]; %#ok<AGROW>
%    PPosSNegPoints = [PPosSNegPoints; n801diffAngles(i) n801threshFull{i}(7)/n801threshFull{i}(2)]; %#ok<AGROW>
% end

for i = 1:length(n557diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n557diffAngles(i) n557threshFull{i}(5)/n557threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n557diffAngles(i) n557threshFull{i}(6)/n557threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n557diffAngles(i) n557threshFull{i}(8)/n557threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n557diffAngles(i) n557threshFull{i}(7)/n557threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n316diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n316diffAngles(i) n316threshFull{i}(5)/n316threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n316diffAngles(i) n316threshFull{i}(6)/n316threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n316diffAngles(i) n316threshFull{i}(8)/n316threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n316diffAngles(i) n316threshFull{i}(7)/n316threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n31diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n31diffAngles(i) n31threshFull{i}(5)/n31threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n31diffAngles(i) n31threshFull{i}(6)/n31threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n31diffAngles(i) n31threshFull{i}(8)/n31threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n31diffAngles(i) n31threshFull{i}(7)/n31threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n801_2diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n801_2diffAngles(i) n801_2threshFull{i}(5)/n801_2threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n801_2diffAngles(i) n801_2threshFull{i}(6)/n801_2threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n801_2diffAngles(i) n801_2threshFull{i}(8)/n801_2threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n801_2diffAngles(i) n801_2threshFull{i}(7)/n801_2threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n558diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n558diffAngles(i) n558threshFull{i}(5)/n558threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n558diffAngles(i) n558threshFull{i}(6)/n558threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n558diffAngles(i) n558threshFull{i}(8)/n558threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n558diffAngles(i) n558threshFull{i}(7)/n558threshFull{i}(2)]; %#ok<AGROW>
end

% least squares linear fits
x = 0:0.001:90;

A = ones(size(PNegSNegPoints));
A(:, 2) = PNegSNegPoints(:,1);
B = squeeze(PNegSNegPoints(:,2));
coeff = A\B;
PNegSNegFit = coeff(1)*ones(size(x)) + coeff(2)*x;

A = ones(size(PNegSPosPoints));
A(:, 2) = PNegSPosPoints(:,1);
B = squeeze(PNegSPosPoints(:,2));
coeff = A\B;
PNegSPosFit = coeff(1)*ones(size(x)) + coeff(2)*x;

A = ones(size(PPosSPosPoints));
A(:, 2) = PPosSPosPoints(:,1);
B = squeeze(PPosSPosPoints(:,2));
coeff = A\B;
PPosSPosFit = coeff(1)*ones(size(x)) + coeff(2)*x;

A = ones(size(PPosSNegPoints));
A(:, 2) = PPosSNegPoints(:,1);
B = squeeze(PPosSNegPoints(:,2));
coeff = A\B;
PPosSNegFit = coeff(1)*ones(size(x)) + coeff(2)*x;

% plotting

figure('position', [100 100 975 450])

axes('position', [0.1 0.15 0.35 0.7], 'fontsize', 18)
hold on

plot(0, 1.9, 's', 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
text(5, 1.9, 'negative on primary electrode', 'fontsize', 16)

plot(0, 1.8, 's', 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
text(5, 1.8, 'positive on primary electrode', 'fontsize', 16)


plot([0 90], [1 1], 'k--')
plot(PNegSNegPoints(:,1), PNegSNegPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
h = plot(x, PNegSNegFit);
set(h, 'color', blue)
plot(PPosSPosPoints(:,1), PPosSPosPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
h = plot(x, PPosSPosFit);
set(h, 'color', rust)

xlabel(['angle between axon', 10, 'and electrode pair (degrees)'], 'fontsize', 20)
ylabel(['normalized threshold'], 'fontsize', 20)
title('same polarity', 'fontsize', 20)
set(gca, 'XLim', [-10 100], 'YLim', [0.4 2], 'XTick', [0 45 90])

axes('position', [0.55 0.15 0.35 0.7], 'fontsize', 18)
hold on
plot([0 90], [1 1], 'k--')
plot(PNegSPosPoints(:,1), PNegSPosPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
h = plot(x, PNegSPosFit);
set(h, 'color', blue)
plot(PPosSNegPoints(:,1), PPosSNegPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
h = plot(x, PPosSNegFit);
set(h, 'color', rust)

title('opposite polarity', 'fontsize', 20)
set(gca, 'XLim', [-10 100], 'YLim', [0.4 2.0], 'XTick', [0 45 90])

% subplot(2,1,2)
% hold on
% plot([0 90], [1 1], 'k--')
% title('ratio of p+, s+ to p+')
% xlabel('angle between axon and electrode pair (degrees)')
% ylabel('threshold ratio')
% set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

% subplot(2,2,4)
% hold on
% plot([0 90], [1 1], 'k--')
% title('ratio of p+, s- to p+')
% set(gca, 'XLim', [-10 100], 'YLim', [0.4 2], 'XTick', [0 45 90])

% subplot(1,2,1); hold on
% plot(PNegSNegPoints(:,1), PNegSNegPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
% h = plot(x, PNegSNegFit);
% set(h, 'color', blue)
% plot(PPosSPosPoints(:,1), PPosSPosPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
% h = plot(x, PPosSPosFit);
% set(h, 'color', rust)
% 
% subplot(1,2,2); hold on
% plot(PNegSPosPoints(:,1), PNegSPosPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
% h = plot(x, PNegSPosFit);
% set(h, 'color', blue)
% plot(PPosSNegPoints(:,1), PPosSNegPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
% h = plot(x, PPosSNegFit);
% set(h, 'color', rust)

% subplot(2,2,3); hold on
% plot(PPosSPosPoints(:,1), PPosSPosPoints(:,2) ,'.', 'MarkerSize', 20)
% plot(x, PPosSPosFit)
% 
% subplot(2,2,4); hold on
% plot(PPosSNegPoints(:,1), PPosSNegPoints(:,2) ,'.', 'MarkerSize', 20)
% plot(x, PPosSNegFit)



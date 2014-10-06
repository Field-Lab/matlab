
errordlg('This script needs to be updated to reflect change in erfFitter to using maximum likelihood')

n227responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/axon_responses_227.csv');
n243responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/axon_responses_243.csv');
n334responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/axon_responses_334.csv');
n558responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/axon_responses_558.csv');
n886responsesFull = csvread('/snle/home/lhruby/Documents/projects/axon_stim_analysis/axon_responses_886.csv');

n334stimAmpsFull = [0.5055 0.5561 0.6066 0.6572 0.7330 0.8088 0.8846 0.9857 1.1043 1.2047 1.3051 1.4055...
    1.6063 1.7067 1.9075 2.1083 2.3091 2.5098 2.8110 3.0118 3.4134 3.7146 4.1161];
n886stimAmpsFull = [0.5055 0.5561 0.6066 0.6572 0.7330 0.8088 0.8846 0.9857 1.1043 1.2047 1.3051 1.4055...
    1.6063 1.7067 1.9075 2.1083 2.3091 2.5098 2.8110 3.0118 3.4134 3.7146];
n243stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043];
n227stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043];
n558stimAmpsFull = [0.1005 0.1131 0.1194 0.1320 0.1445 0.1634 0.1759 0.1948 0.2136 0.2388 0.2576 0.2780...
    0.3033 0.3539 0.3791 0.4297 0.4550 0.5055 0.5561 0.6066 0.6824 0.7330 0.8088 0.8846 0.9857 1.1043...
    1.2047 1.3051 1.4055 1.6063 1.7067 1.9075 2.1083];

n227patternNos = n227responsesFull(:,1);
n243patternNos = n243responsesFull(:,1);
n334patternNos = n334responsesFull(:,1);
n558patternNos = n558responsesFull(:,1);
n886patternNos = n886responsesFull(:,1);
n227responses = cell(length(n227patternNos), 1);
n243responses = cell(length(n243patternNos), 1);
n334responses = cell(length(n334patternNos), 1);
n558responses = cell(length(n558patternNos), 1);
n886responses = cell(length(n886patternNos), 1);
n227stimAmps = cell(length(n227patternNos), 1);
n243stimAmps = cell(length(n243patternNos), 1);
n334stimAmps = cell(length(n334patternNos), 1);
n558stimAmps = cell(length(n558patternNos), 1);
n886stimAmps = cell(length(n886patternNos), 1);
n227params = cell(length(n227patternNos), 1);
n243params = cell(length(n243patternNos), 1);
n334params = cell(length(n334patternNos), 1);
n558params = cell(length(n558patternNos), 1);
n886params = cell(length(n886patternNos), 1);
n227projections = cell(length(n227patternNos), 1);
n243projections = cell(length(n243patternNos), 1);
n334projections = cell(length(n334patternNos), 1);
n558projections = cell(length(n558patternNos), 1);
n886projections = cell(length(n886patternNos), 1);
n227thresh = zeros(length(n227patternNos), 1);
n243thresh = zeros(length(n243patternNos), 1);
n334thresh = zeros(length(n334patternNos), 1);
n558thresh = zeros(length(n558patternNos), 1);
n886thresh = zeros(length(n886patternNos), 1);

for i = 1:length(n227patternNos)
    n227responses{i} = n227responsesFull(i, 2:end);
    n227stimAmps{i} = n227stimAmpsFull;
    for j = length(n227responses{i}):-1:1
        if n227responses{i}(j) > 1.1
            n227responses{i}(j) = [];
            n227stimAmps{i}(j) = [];
        end
    end
    [n227params{i} n227projections{i}] = erfFitter([n227stimAmps{i}; n227responses{i}], 1, 1);
    n227thresh(i) = -n227params{i}(2)/n227params{i}(1);
end

for i = 1:length(n243patternNos)
    n243responses{i} = n243responsesFull(i, 2:end);
    n243stimAmps{i} = n243stimAmpsFull;
    for j = length(n243responses{i}):-1:1
        if n243responses{i}(j) > 1.1
            n243responses{i}(j) = [];
            n243stimAmps{i}(j) = [];
        end
    end
    [n243params{i} n243projections{i}] = erfFitter([n243stimAmps{i}; n243responses{i}], 1, 1);
    n243thresh(i) = -n243params{i}(2)/n243params{i}(1);
end


for i = 1:length(n334patternNos)
    n334responses{i} = n334responsesFull(i, 2:end);
    n334stimAmps{i} = n334stimAmpsFull;
    for j = length(n334responses{i}):-1:1
        if n334responses{i}(j) > 1.1
            n334responses{i}(j) = [];
            n334stimAmps{i}(j) = [];
        end
    end
    [n334params{i} n334projections{i}] = erfFitter([n334stimAmps{i}; n334responses{i}], 1, 1);
    n334thresh(i) = -n334params{i}(2)/n334params{i}(1);
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


for i = 1:length(n886patternNos)
    n886responses{i} = n886responsesFull(i, 2:end);
    n886stimAmps{i} = n886stimAmpsFull;
    for j = length(n886responses{i}):-1:1
        if n886responses{i}(j) > 1.1
            n886responses{i}(j) = [];
            n886stimAmps{i}(j) = [];
        end
    end
    [n886params{i} n886projections{i}] = erfFitter([n886stimAmps{i}; n886responses{i}], 1, 1);
    n886thresh(i) = -n886params{i}(2)/n886params{i}(1);
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

%% plotting: neuron 227 (on parasol)
movies = []; %for now, a required argument to plot, but not used

dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 227;
responses = n227responses;
erfFitDetails.stimAmps    = n227stimAmps;
erfFitDetails.projections = n227projections;
erfFitDetails.thresholds  = n227thresh;
erfFitDetails.patternNos  = n227patternNos;

nPairs = 2;
n227threshFull = cell(nPairs, 1);
n227pairIDs = cell(nPairs, 1);

% electrode pair 35, 37 (primary = 35)
pElec = 35;
sElec = 37;
patterns = [49 50 57 58 13:16];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n227threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end


% % electrode pair 35, 39 (primary = 35)
% pElec = 35;
% sElec = 39;
% patterns = [49 50 61 62 21:24];
% axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)
% 
% for i = 1:8
%     patternIndex = find(erfFitDetails.patternNos == patterns(i));
%     n227threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
% end
% n227pairIDs{2} = [pElec sElec];

% electrode pair 37, 39 (primary = 37)
pElec = 37;
sElec = 39;
patterns = [57 58 61 62 41:44];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n227threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
end
n227pairIDs{2} = [pElec sElec];

%not enough of the resopnse curve is measured to get threshold, so can only use lower bound
n227threshFull{2}(6) = 1.1043;
n227threshFull{2}(7) = 1.1043;


%% plotting: neuron 243 (off midget)
dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data007';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data006/data006.ei';
neuronID = 243;
responses = n243responses;
erfFitDetails.stimAmps    = n243stimAmps;
erfFitDetails.projections = n243projections;
erfFitDetails.thresholds  = n243thresh;
erfFitDetails.patternNos  = n243patternNos;

nPairs = 1;
n243threshFull = cell(nPairs, 1);
n243pairIDs = cell(nPairs, 1);

% electrode pair 35, 33 (primary = 35)
pElec = 35;
sElec = 33;
patterns = [49:50 53:54 5:8];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n243threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end
n243pairIDs{1} = [pElec sElec];

%% plotting: neuron 334 (off parasol)
dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data011';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data012/data012.ei';
neuronID = 334;
responses = n334responses;
erfFitDetails.stimAmps    = n334stimAmps;
erfFitDetails.projections = n334projections;
erfFitDetails.thresholds  = n334thresh;
erfFitDetails.patternNos  = n334patternNos;

nPairs = 1;
n334threshFull = cell(nPairs, 1);
n334pairIDs = cell(nPairs, 1);

% % electrode pair 56, 60 (primary = 60) % no responses to 56 alone and little effect on 60
% pElec = 60;
% sElec = 56;
% patterns = [121 122 111 112 79 81 80 82];
% axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

% electrode pair 60, 61 (primary = 61)
pElec = 61;
sElec = 60;
patterns = [123 124 121 122 107 109 108 110];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n334threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end
n334pairIDs{1} = [pElec sElec];

%% plotting: neuron 558 (on parasol)
dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data010';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data009/data009.ei';
neuronID = 558;
responses = n558responses;
erfFitDetails.stimAmps    = n558stimAmps;
erfFitDetails.projections = n558projections;
erfFitDetails.thresholds  = n558thresh;
erfFitDetails.patternNos  = n558patternNos;

nPairs = 7;
n558threshFull = cell(nPairs, 1);
n558pairIDs = cell(nPairs, 1);

% electrode pair 3, 1 (primary = 3)
pElec = 3;
sElec = 1;
patterns = [49 50 51 52 1:4];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{1} = [pElec sElec];

% electrode pair 3, 2 (primary = 2)
pElec = 2;
sElec = 3;
patterns = [53 54 49 50 5 7 6 8];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{2}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{2} = [pElec sElec];

% electrode pair 2, 5 (primary = 2)
pElec = 2;
sElec = 5;
patterns = [53 54 55 56 33:36];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{3}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{3} = [pElec sElec];

% electrode pair 3, 5 (primary = 5)
pElec = 3;
sElec = 5;
patterns = [49 50 55 56 9:12];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{4}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{4} = [pElec sElec];

% electrode pair 3, 6 (primary = 3)
pElec = 3;
sElec = 6;
patterns = [49 50 57 58 13:16];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{5}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{5} = [pElec sElec];

% electrode pair 3, 7 (primary = 3)
pElec = 3;
sElec = 7;
patterns = [49 50 59 60 17:20];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{6}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{6} = [pElec sElec];

% electrode pair 3, 64 (primary = 3)
pElec = 3;
sElec = 64;
patterns = [49 50 61 62 21:24];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)

for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n558threshFull{7}(i) = erfFitDetails.thresholds(patternIndex);
end
n558pairIDs{7} = [pElec sElec];

%% plotting: neuron 886 (on parasol)
dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data011';
pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data012/data012.ei';
neuronID = 886;
responses = n886responses;
erfFitDetails.stimAmps    = n886stimAmps;
erfFitDetails.projections = n886projections;
erfFitDetails.thresholds  = n886thresh;
erfFitDetails.patternNos  = n886patternNos;

nPairs = 1;
n886threshFull = cell(nPairs, 1);
n886pairIDs = cell(nPairs, 1);

% electrode pair 60, 61 (primary = 60)
pElec = 60;
sElec = 61;
patterns = [121 122 123 124 107 108 109 110];
axonStimFigurePlotter(responses, pElec, sElec, patterns, movies, dataPath, pathToEi, neuronID, erfFitDetails)


for i = 1:8
    patternIndex = find(erfFitDetails.patternNos == patterns(i));
    n886threshFull{1}(i) = erfFitDetails.thresholds(patternIndex);
end
n886pairIDs{1} = [pElec sElec];

%% defining electrode angles (vertex = primary, angle drawn from positive x axis as plotted)
n227pairAngles = [210 90];
n243pairAngles = 30;
n334pairAngles = 270;
n558pairAngles = [210 270 330 30 270 330 150];
n886pairAngles = 90;

%% defining angle of axon (vertex = soma, angle draw from positive x axis as plotted)
n227axonAngle = 186.6;
n243axonAngle = 188.8;
n334axonAngle = 109.3;
n558axonAngle = 31.3;
n886axonAngle = 111.2;

n227diffAngles = abs(n227pairAngles - ones(size(n227pairAngles))*n227axonAngle);
n243diffAngles = abs(n243pairAngles - ones(size(n243pairAngles))*n243axonAngle);
n334diffAngles = abs(n334pairAngles - ones(size(n334pairAngles))*n334axonAngle);
n558diffAngles = abs(n558pairAngles - ones(size(n558pairAngles))*n558axonAngle);
n886diffAngles = abs(n886pairAngles - ones(size(n886pairAngles))*n886axonAngle);

%% converts angles to equivalent within (1, 90)

for i = 1:length(n227diffAngles)
    if n227diffAngles(i) > 180;
        n227diffAngles(i) = 360 - n227diffAngles(i);
    end
    if n227diffAngles(i) > 90
        n227diffAngles(i) = 180 - n227diffAngles(i);
    end
end

for i = 1:length(n243diffAngles)
    if n243diffAngles(i) > 180;
        n243diffAngles(i) = 360 - n243diffAngles(i);
    end
    if n243diffAngles(i) > 90
        n243diffAngles(i) = 180 - n243diffAngles(i);
    end
end

for i = 1:length(n334diffAngles)
    if n334diffAngles(i) > 180;
        n334diffAngles(i) = 360 - n334diffAngles(i);
    end
    if n334diffAngles(i) > 90
        n334diffAngles(i) = 180 - n334diffAngles(i);
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

for i = 1:length(n886diffAngles)
    if n886diffAngles(i) > 180;
        n886diffAngles(i) = 360 - n886diffAngles(i);
    end
    if n886diffAngles(i) > 90
        n886diffAngles(i) = 180 - n886diffAngles(i);
    end
end



%% summary plot

PNegSNegPoints = [];
PNegSPosPoints = [];
PPosSNegPoints = [];
PPosSPosPoints = [];

for i = 1:length(n227diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n227diffAngles(i) n227threshFull{i}(5)/n227threshFull{i}(1)]; %#ok<AGROW>
   if i == 1
       PNegSPosPoints = [PNegSPosPoints; n227diffAngles(i) n227threshFull{i}(6)/n227threshFull{i}(1)]; %#ok<AGROW>
   end   
   PPosSPosPoints = [PPosSPosPoints; n227diffAngles(i) n227threshFull{i}(8)/n227threshFull{i}(2)]; %#ok<AGROW>   
   if i == 1
       PPosSNegPoints = [PPosSNegPoints; n227diffAngles(i) n227threshFull{i}(7)/n227threshFull{i}(2)]; %#ok<AGROW>
   end
end

for i = 1:length(n243diffAngles)

   PNegSNegPoints = [PNegSNegPoints; n243diffAngles(i) n243threshFull{i}(5)/n243threshFull{i}(1)]; %#ok<AGROW>   
   PNegSPosPoints = [PNegSPosPoints; n243diffAngles(i) n243threshFull{i}(6)/n243threshFull{i}(1)]; %#ok<AGROW>   
   PPosSPosPoints = [PPosSPosPoints; n243diffAngles(i) n243threshFull{i}(8)/n243threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n243diffAngles(i) n243threshFull{i}(7)/n243threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n334diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n334diffAngles(i) n334threshFull{i}(5)/n334threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n334diffAngles(i) n334threshFull{i}(6)/n334threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n334diffAngles(i) n334threshFull{i}(8)/n334threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n334diffAngles(i) n334threshFull{i}(7)/n334threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n558diffAngles)

   PNegSNegPoints = [PNegSNegPoints; n558diffAngles(i) n558threshFull{i}(5)/n558threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n558diffAngles(i) n558threshFull{i}(6)/n558threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n558diffAngles(i) n558threshFull{i}(8)/n558threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n558diffAngles(i) n558threshFull{i}(7)/n558threshFull{i}(2)]; %#ok<AGROW>
end

for i = 1:length(n886diffAngles)
   PNegSNegPoints = [PNegSNegPoints; n886diffAngles(i) n886threshFull{i}(5)/n886threshFull{i}(1)]; %#ok<AGROW>
   PNegSPosPoints = [PNegSPosPoints; n886diffAngles(i) n886threshFull{i}(6)/n886threshFull{i}(1)]; %#ok<AGROW>
   PPosSPosPoints = [PPosSPosPoints; n886diffAngles(i) n886threshFull{i}(8)/n886threshFull{i}(2)]; %#ok<AGROW>
   PPosSNegPoints = [PPosSNegPoints; n886diffAngles(i) n886threshFull{i}(7)/n886threshFull{i}(2)]; %#ok<AGROW>
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

%% plotting (new version)

blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon


% plotting
figure('position', [100 100 975 450])

axes('position', [0.1 0.15 0.35 0.7], 'fontsize', 18)
hold on

plot(0, 1.51, 's', 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
text(5, 1.51, '   on primary electrode', 'fontsize', 16)

plot(0, 1.42, 's', 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
text(5, 1.42, '   on primary electrode', 'fontsize', 16)


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
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

axes('position', [0.55 0.15 0.35 0.7], 'fontsize', 18)
hold on
plot([0 90], [1 1], 'k--')
plot(PNegSPosPoints(:,1), PNegSPosPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue)
h = plot(x, PNegSPosFit);
set(h, 'color', blue)
plot(PPosSNegPoints(:,1), PPosSNegPoints(:,2) ,'.', 'MarkerSize', 20, 'MarkerFaceColor', rust, 'MarkerEdgeColor', rust)
h = plot(x, PPosSNegFit);
set(h, 'color', rust)

plot(n227diffAngles(2), n227threshFull{2}(6)/n227threshFull{2}(1),'+', 'MarkerSize', 10, 'MarkerEdgeColor', blue)
plot(n227diffAngles(2), n227threshFull{2}(7)/n227threshFull{2}(2),'+', 'MarkerSize', 10, 'MarkerEdgeColor', rust)


title('opposite polarity', 'fontsize', 20)
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])



%% plotting (old version)

figure
subplot(2,2,1)
hold on
plot([0 90], [1 1], 'k--')
title('ratio of p-, s- to p-')
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

subplot(2,2,2)
hold on
plot([0 90], [1 1], 'k--')
title('ratio of p-, s+ to p-')
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

subplot(2,2,3)
hold on
plot([0 90], [1 1], 'k--')
title('ratio of p+, s+ to p+')
xlabel('angle between axon and electrode pair (degrees)')
ylabel('threshold ratio')
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

subplot(2,2,4)
hold on
plot([0 90], [1 1], 'k--')
title('ratio of p+, s- to p+')
set(gca, 'XLim', [-10 100], 'YLim', [0.4 1.6], 'XTick', [0 45 90])

subplot(2,2,1); hold on
plot(PNegSNegPoints(:,1), PNegSNegPoints(:,2) ,'.', 'MarkerSize', 20)
plot(x, PNegSNegFit)

subplot(2,2,2); hold on
plot(PNegSPosPoints(:,1), PNegSPosPoints(:,2) ,'.', 'MarkerSize', 20)
plot(n227diffAngles(2), n227threshFull{2}(6)/n227threshFull{2}(1),'+', 'MarkerSize', 10)
plot(x, PNegSPosFit)

subplot(2,2,3); hold on
plot(PPosSPosPoints(:,1), PPosSPosPoints(:,2) ,'.', 'MarkerSize', 20)
plot(x, PPosSPosFit)

subplot(2,2,4); hold on
plot(PPosSNegPoints(:,1), PPosSNegPoints(:,2) ,'.', 'MarkerSize', 20)
plot(n227diffAngles(2), n227threshFull{2}(7)/n227threshFull{2}(2),'+', 'MarkerSize', 10)
plot(x, PPosSNegFit)


%% collecting the single-electrode stimulation patterns


singleElecAxonResponses{1} = n227responses{37};
singleElecAxonProjections{1} = n227projections{37};
singleElecAxonThresh(1) = n227thresh(37);
singleElecAxonStimAmps{1} = n227stimAmps{37};

singleElecAxonResponses{2} = n243responses{5};
singleElecAxonProjections{2} = n243projections{5};
singleElecAxonThresh(2) = n243thresh(5);
singleElecAxonStimAmps{2} = n243stimAmps{5};

singleElecAxonResponses{3} = n334responses{7};
singleElecAxonProjections{3} = n334projections{7};
singleElecAxonThresh(3) = n334thresh(7);
singleElecAxonStimAmps{3} = n334stimAmps{7};

singleElecAxonResponses{4} = n558responses{41};
singleElecAxonProjections{4} = n558projections{41};
singleElecAxonThresh(4) = n558thresh(41);
singleElecAxonStimAmps{4} = n558stimAmps{41};

singleElecAxonResponses{5} = n886responses{5};
singleElecAxonProjections{5} = n886projections{5};
singleElecAxonThresh(5) = n886thresh(5);
singleElecAxonStimAmps{5} = n886stimAmps{5};











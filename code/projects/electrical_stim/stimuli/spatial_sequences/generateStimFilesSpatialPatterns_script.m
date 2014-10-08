% child functions reworked and tested November 2010

%2009-09-03-1 data003-data008
% electrodes = [7 17 23 34 44]; %4 electrodes that have 4 different responding neurons that each reach 100% response
% fullAmps = [0.3 0.9 1.5 1 1.3]; %amplitudes that give 100% response in corresponding target neurons
% threshAmps = [0.08 0.5 0.5 0.6 0.7]; %amplitudes that give 50% response in corresponding target neurons

%2010-10-18-3 data023
%electrodes = [17 21 47 36];
%fullAmps = [0.35 1.1 0.8 0.5];

% electrodes = [1 2 3 4];
% fullAmps = [1.2 1.4 1.6 1.8];
% threshAmps = [0.6 0.7 0.8 0.9];

% electrodes = [14  39   50    63];
% threshAmps = [1   0.9  1.1   0.7];
% fullAmps =   [1.6 1.3  1.65  1];

% electrodes = [60   14   21    26];
% threshAmps = [0.6  0.4  0.5   0.7];
% fullAmps =   [0.9  0.7  0.7   0.9];

% electrodes = [23   31    33    40];
% threshAmps = [1    0.4   0.7   1.1];
% fullAmps =   [1.5  0.6   1.1   1.7];

%electrodes = [23   31    33    40];
%threshAmps = [1    0.4   0.7   1.1];
%fullAmps =   [1.5  0.6   1.1   1.7];

% electrodes = [15   5    11    36    63];
% threshAmps = [1    1.2  1     0.9   0.65];
% fullAmps =   [1.3  1.8  1.4   1.15  0.9];

% electrodes = [4    16   24    28    54   61];
% threshAmps = [0.7  0.6  0.35  0.25  0.8  0.37];
% fullAmps =   [1    0.8  0.6   0.35  1.15 0.6];

electrodes = [18    22    49    38];
threshAmps = [1.1   2.2   2.4   1.7];
fullAmps   = [1.8   3     3.5   2.4];

periodInMs = 50; %time between application of difference sequences -- make sure there is enough time to play entire sequence
offsets = [1 2 5]; %delays between subsequent pulses in sequences (in ms) -- must be 3 values
nOrders = 4; %number of random orders of pulses to generate


generateStimFilesSpatialPatterns(electrodes, fullAmps, threshAmps, offsets, nOrders, periodInMs)
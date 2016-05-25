%% Example 4: Use stimulating electrode information (still not advisable, very unstable)
%Gonzalo Mena, 3/16



%% Comments

% This example does spike sorting for pattern 116 (equal to stimulating
% electrode), and neurons 4058 4656 3457.
% Preparation: 2012-09-24-3
% Pattern is in folder data004 (no need to specify, but if not, one is at
% risk that the algorithm will look for it in other folders.
% alternatively, one can specify a subset of folder indexes to look for
% (see below).


%%
%Initialization: Do this only once, a 'profile' of the artifact will be
%taken based on a single pattern. One can specify this pattern/folder index (e.g. 4
%for data004), only pattern of none of them. (see below)
%Initialization should take between 2 and 10 minutes, and again: need to do
%it once for each preparation. Then, only care about variable params.

patternNoInitial=200;  %Pattern chosen to do initialization. Can leave empty and will use the first it finds
                       %Ideally, chose one with little axonal bundle, and
                       %as closest to the center of the array as possible,
                       %to avoid undesided 'aliasing' effects.
                       
IndFolderInitial={'data003','data004','data005','data006'};  %look for the pattern in data003,data004,data005,data006. If empty, will look at all folderes
pathToPreparationInitial='/Volumes/Analysis/2012-09-24-3/';

params=InitializeArray(pathToPreparationInitial,patternNoInitial);



%% Change some default values, 

params.global.tarray=[0 [7:30]]; %(times to look for spikes);
params.global.thresEI=35; %(for spike sorting, only consider electrodes with strong enough signal)
params.global.Tmax=40; %use first 2ms of recordings (recall, sampling rate 20KhZ)
params.global.nTrials=50; %maximum number of trials for the same stimulus (this is to create appropriate matrices, unbalances
                           %are solved by filling with NaNs (see TracesAll)
params.bundle.findBundle=1;
params.bundle.cutBundle=1;  %Terminate spike sorting  just before the onset of bundle
params.bundle.useBundleAlg=0; %dont change the variance model after the bundle
params.global.useStimElec=1;
params.global.sortData=1;

%%Set values for spike sorting, and do spike sorting
neuronIds=[4058 4656 3457];
patternNo=116;
IndFolder={'data003','data004','data005','data006'};  %look for the pattern in data003,data004,data005,data006. If empty, will look at all folderes
                      %it is optional, but if not stated, it may use
                      %information of an undesired folder.
                      
pathToEi=['/Volumes/Analysis/data000/data000.ei'];
pathToPreparation='/Volumes/Analysis/2012-09-24-3/';
tic
[Output]=DoSpikeSortingLargeScale(pathToPreparation,pathToEi,patternNo,neuronIds,params,IndFolder);
toc
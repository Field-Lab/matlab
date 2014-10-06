% Code written by Nishal

% To test spike sorting on data in lab
%% Load paths
startup
%% Load raw File
   % open raw data file
   
   
%   rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile('/Volumes/Data/2005-04-26-1/data006/data006000.bin');
% 
%     % gets data from sample 0 to 20000 (first second) across all electrodes
%     samplingRate = 20000;
%     time_length=72;
%     data = rawFile.getData(0, time_length*samplingRate);
% 
%     % Stimulus TTLs (most times every 100 frames) on java electrode position 0 (1 in matlab)
%   
% 
%     rawFile.close();


% open raw data file 
    rawFile = edu.ucsc.neurobiology.vision.io.RawDataWrapper('/Volumes/Data/2005-04-26-1/data006');

    samplingRate = 20000;

    %set up arguments
    sample = 1000; %sample to start reading from
    nSamples = 72*samplingRate; % get one second, starting from sample
    nElectrodes = 512; % ALL electrodes in the dataset

    %get the raw data
    data = rawFile.getData(sample, nSamples, nElectrodes);

    % Stimulus TTLs (most times every 100 frames) on java electrode position 0 (1 in matlab)
   
    rawFile.close();
    data=data';
    data=double(data);
    electrode_list=[31,23,39,35,27,36,28];
    data = data(electrode_list,:);
    dt = 1/samplingRate;
    
   clear rawFile timelength
   
   save('../CBPSpikesortDemoPackage/example_data/2005-04-26-1_nps.mat');
   
   %% 
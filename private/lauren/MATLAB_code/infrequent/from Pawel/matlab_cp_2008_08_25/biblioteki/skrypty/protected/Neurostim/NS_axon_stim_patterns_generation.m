TimeShiftInMs=0;
DelayInMs=30;

NumberOfSamples=40000;

CenterElectrodes=[39];
nLayers = 1; %use 1 for 7-electrode cluster and 2 for 42-electrode cluster
NumberOfClusters=length(CenterElectrodes);

[electrodes Array] = generateAxonStimPatterns(CenterElectrodes,nLayers);

[PatternsOut,Times,MovieChunk]=NS_MovieChunksForExperimentAxon(NumberOfClusters,Array,TimeShiftInMs,DelayInMs,NumberOfSamples);

MovieChunksFile=[1 MovieChunk]; %only one movie

%break;
cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('axon_electrodes','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('axon_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('axon_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
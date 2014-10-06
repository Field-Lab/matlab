CenterElectrodes=[60];
TimeShiftInMs=0;
DelayInMs=30;
NumberOfSamples=50000;

nLayers=2;
[electrodes Array] = generateAxonStimPatterns(CenterElectrodes,1,nLayers);
[PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2elClusterOf19(NumberOfClusters,TimeShiftInMs,DelayInMs,NumberOfSamples);

cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('2el_electrodes_19','wb','ieee-le.l64')
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('2el_patterns_19','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('2el_movie_19','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
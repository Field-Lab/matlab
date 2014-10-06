CenterElectrodes=[3 19 55];
TimeShiftInMs=0;
DelayInMs=30;
NumberOfSamples=30000;

[electrodes Array] = generateClusterStimPatterns2ElScan(CenterElectrodes);

[PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2el(length(CenterElectrodes),TimeShiftInMs,DelayInMs,NumberOfSamples);
% break
disp('?')
% cd C:\pawel\pliki\nauka\matlab\; 
cd ~/Desktop/StimFiles
fid = fopen('2el_electrodes','wb')
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('2el_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('2el_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
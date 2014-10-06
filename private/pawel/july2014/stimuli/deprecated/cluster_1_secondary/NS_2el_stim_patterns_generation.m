CenterElectrodes=[3 18 50];
TimeShiftInMs=0;
DelayInMs=30;
NumberOfSamples=30000;

[electrodes Array] = generateClusterStimPatterns2ElScan(CenterElectrodes);


%%% NS_MovieChunksForExperiment2el is now equivalent to modified version (edited by Pawel)
[PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2el(length(CenterElectrodes),TimeShiftInMs,DelayInMs,NumberOfSamples);

keyboard
%break;
cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('2el_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('2el_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double'); %was previously commented out for some reason...
fclose(fid);

fid = fopen('2el_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
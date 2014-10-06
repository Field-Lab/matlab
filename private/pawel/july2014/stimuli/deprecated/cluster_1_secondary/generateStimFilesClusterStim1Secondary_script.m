clear all

%%%%%%%%%%%
%This script and associated function have not been fully tested!!!!!
%%%%%%%%%%%

centerElectrodes=[8 20];
TimeShiftInMs=0;
DelayInMs=30;
%NumberOfSamples=30000;
relAmps = [-0.2 -0.1 0.1 0.2];
maxDistance = 1.8;

%[electrodes Array] = generateClusterStimPatterns2ElScan(CenterElectrodes);

[electrodes Array clusterIDs] = generatePatternClusterStim1SecondaryWrapper(centerElectrodes, relAmps, maxDistance);

%[PatternOrder,MovieChunksFile]=NS_MovieChunksForExperiment2el(length(CenterElectrodes),TimeShiftInMs,DelayInMs,NumberOfSamples);
[PatternOrder, MovieChunksFile] = generateMovieClusterStim1Secondary(length(centerElectrodes), clusterIDs, TimeShiftInMs, DelayInMs);

fid = fopen('cluster1sec_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('cluster1sec_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double'); %was previously commented out for some reason...
fclose(fid);

fid = fopen('cluster1sec_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
% generates STIM64 stimulus files for stimulating with 3-5 electrodes simultaneously and at small
% time offsets (1,5,10 ms)


electrodes=[3 18 38 22]; %4 electrodes that have 4 different responding neurons that each reach 100% response
stimAmpsFullResponse = [0.5 0.4 0.9 1]; %amplitudes that give 100% response in corresponding target neurons
stimAmpsThresh = [0.3 0.2 0.6 0.8]; %threshold amplitudes for corresponding target neurons

TimeShiftInMs=0;
DelayInMs=100; %time between application of patterns (in ms)
NumberOfSamples=24600; %length of each movie chunk (samples)
% refresh rate in labview must be the same or larger than this number!!! (limit is 40000 = 2 s)

[electrodes Array] = generateSpatialPatternStimPatterns(electrodes, stimAmpsFullResponse, stimAmpsThresh);

[PatternOrder,MovieChunksFile]=NS_MovieChunksForExperimentSpatialPatterns(length(electrodes), TimeShiftInMs, DelayInMs, NumberOfSamples);

fid = fopen('spatial_electrodes','wb');
fwrite(fid, electrodes, 'int32');
fclose(fid);

fid = fopen('spatial_patterns','wb','ieee-le.l64');
fwrite(fid, Array, 'double'); %was previously commented out for some reason...
fclose(fid);

fid = fopen('spatial_movie','wb','ieee-le.l64');
fwrite(fid, MovieChunksFile, 'int32');
fclose(fid);
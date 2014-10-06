TimeShiftInMs=0;
DelayInMs=0;
NumberOfSamples=20000;

PatternNumber=6;
Periods=[200 100 50 25 12.5 6.25]*20;
MovieChunksFile=NS_MovieChunksForExperimentFrequencyScan(PatternNumber,Periods,DelayInMs,NumberOfSamples);

fid = fopen('FreqScan_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
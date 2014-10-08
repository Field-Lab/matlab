TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000; %trigger interval = 0.5 s
electrodes=[1:512];
Array=eye(512,512); 

Patterns=NS512_OptimalElectrodeSequence();

Times=[TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];

% cd /Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_files_test; 
% cd /Users/grosberg/Desktop/  


AllPatterns=[];

Chunks=[];
for i=1:8 % 8 Movie ChunksA
    Start=(i-1)*64;
    PatternsForMovie=Patterns(Start+1:Start+64);
    mChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
%     Chunks=[Chunks Chunk];
    AllPatterns=[AllPatterns mChunk];
end

MovieChunksFile=[8 AllPatterns];

keyboard;

fid=fopen('512el_mv','wb');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);

fid = fopen('512el_el','wb');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512el_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

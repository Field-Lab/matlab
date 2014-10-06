TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

Patterns=NS512_OptimalElectrodeSequence();

Times=[TimeShift:InterPulseLatency:TimeShift+31*InterPulseLatency];

Chunks=[];
for i=1:16
    PatternsForMovie=Patterns((i-1)*32+1:i*32)
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end
%break;
MovieChunksFile=[8 Chunks]; %only one movie
%break;
cd C:\home\pawel\2010\stim_files; 
fid = fopen('512_16mv_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512_16mv_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('512_16mv_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
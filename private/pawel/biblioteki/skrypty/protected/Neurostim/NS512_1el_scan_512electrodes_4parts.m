TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

Patterns=NS512_OptimalElectrodeSequence();

Times=[TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];

cd C:\home\pawel\2012\stim_files; 
           
fid = fopen('128el_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('128el_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

Chunks=[];
for i=1:4
    for j=1:2
        Start=(i-1)*128+(j-1)*64;
        PatternsForMovie=Patterns(Start+1:Start+64)
        Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
        Chunks=[Chunks Chunk];
    end
    MovieChunksFile=[2 Chunks];
    name=['128el_mv' num2str(i)];
    fid=fopen(name,'wb');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid); 
end
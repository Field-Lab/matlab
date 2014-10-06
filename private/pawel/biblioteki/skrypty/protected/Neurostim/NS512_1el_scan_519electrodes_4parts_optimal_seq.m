TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

Patterns=NS512_519OptimalOrder
Patterns1=Patterns
Electrodes=NS512_519InverseMess(Patterns);
Patterns=Electrodes;

Times=[TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];

cd /Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_files_test; 
           
AllPatterns=[];

for i=1:4
    Chunks=[];
    StartIndex=(i-1)*128;
    PatternsForSet=zeros(1,128);
    PatternsForSet(1:32)=Patterns(StartIndex+1:4:StartIndex+128);
    PatternsForSet(33:64)=Patterns(StartIndex+2:4:StartIndex+128);
    PatternsForSet(65:96)=Patterns(StartIndex+3:4:StartIndex+128);
    PatternsForSet(97:128)=Patterns(StartIndex+4:4:StartIndex+128);
    
    for j=1:2
        Start=(j-1)*64;
        PatternsForMovie=PatternsForSet(Start+1:Start+64);
        Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
        Chunks=[Chunks Chunk];
        AllPatterns=[AllPatterns PatternsForMovie]; % just to check
    end
    MovieChunksFile=[2 Chunks];
    name=['128el519_mv' num2str(i)];
    fid=fopen(name,'wb');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid); 
end

fid = fopen('128el519_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('128el519_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid); 
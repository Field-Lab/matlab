TimeShiftInMs=15;
InterPulseLatencyInMs=30;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];

Chunks=[];

Rows=[1:16];
WhichColumns=[1:17];
Array1=zeros(length(Rows),512);
ble=NS512_PatternsForRowStimulation(Rows,WhichColumns);
for i=1:length(Rows)
    Array1(i,ble(i,:))=1;
end
Patterns=1:length(Rows);
Times=[TimeShift:InterPulseLatency:TimeShift+length(Rows)*InterPulseLatency];
Chunk1=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);

Rows2=[1:16];
WhichColumns=[18:32];
Array2=zeros(length(Rows2),512);
ble=NS512_PatternsForRowStimulation(Rows2,WhichColumns);
for i=1:length(Rows2)
    Array2(i,ble(i,:))=1;
end
Patterns2=17:length(Rows)+length(Rows2);
Times=[TimeShift:InterPulseLatency:TimeShift+length(Rows2)*InterPulseLatency];
Chunk2=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns2);

Array=[Array1' Array2']';
MovieChunksFile=[2 Chunk1 Chunk2]; %only one movie

cd C:\home\pawel\2010\stim_files; 

fid = fopen('512_Rows_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512_Rows_pt','wb','ieee-le.l64')
fwrite(fid,Array','double');
fclose(fid);

fid = fopen('512_Rows_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
TimeShiftInMs=15;
InterPulseLatencyInMs=30;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];

Chunks=[];

Columns=[1:16];
WhichRows=[1:16];
Array1=zeros(length(Columns),512);
ble=NS512_PatternsForColumnStimulation(Columns,WhichRows);
for i=1:length(Columns)
    Array1(i,ble(i,:))=1;
end
Patterns=1:length(Columns);
Times=[TimeShift:InterPulseLatency:TimeShift+length(Columns)*InterPulseLatency];
Chunk1=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);

Columns=[32:-1:17];
WhichRows=[1:16];
Array2=zeros(length(Columns),512);
ble=NS512_PatternsForColumnStimulation(Columns,WhichRows);
for i=1:length(Columns)
    Array2(i,ble(i,:))=1;
end
Patterns=1:length(Columns);
Times=[TimeShift:InterPulseLatency:TimeShift+length(Columns)*InterPulseLatency];
Chunk2=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns+16);

Array=[Array1' Array2']';
MovieChunksFile=[2 Chunk1 Chunk2]; %only one movie

%break;
cd C:\home\pawel\2010\stim_files; 
fid = fopen('512_Culumns_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512_Columns_pt','wb','ieee-le.l64')
fwrite(fid,Array','double');
fclose(fid);

fid = fopen('512_Columns_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
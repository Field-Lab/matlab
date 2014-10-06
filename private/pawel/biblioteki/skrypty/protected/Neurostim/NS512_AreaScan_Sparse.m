TimeShiftInMs=15;
InterPulseLatencyInMs=15;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];

PatternsSequence=NS512_OptimalElectrodeSequence();

RowsIndexes=[1:2:16];
ColumnsIndexes=[17:2:32];
Electrodes=NS512_PatternsForRectangularAreaStimulation(RowsIndexes,ColumnsIndexes);

indexes=[];
for i=1:length(Electrodes)       
    electrode=Electrodes(i);
    index=find(PatternsSequence==electrode);
    indexes=[indexes index];
    a=sort(indexes);
    Patterns=PatternsSequence(a)
end
    
%break;
Array=eye(512,512); 

Times=[TimeShift:InterPulseLatency:TimeShift+31*InterPulseLatency];

Chunks=[];
for i=1:2
    PatternsForMovie=Patterns((i-1)*32+1:i*32);
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end
    
%Patterns=ones(1,1);

MovieChunksFile=[2 Chunks]; %only one movie

%break;
cd F:\StimFiles; 

fid = fopen('512_1el_area_sparse_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512_1el_area_sparse_pt','wb','ieee-le.l64')
fwrite(fid,Array','double');
fclose(fid);

fid = fopen('512_1el_area_sparse_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
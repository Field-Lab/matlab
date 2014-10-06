TimeShiftInMs=30;
DelayInMs=30;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=zeros(512,64); 
Array=eye(512,512); 
Times=[];

Patterns=NS512_OptimalElectrodeSequence();
Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);
StartColumn=1;
StartRow=1;

Columns=[StartColumn:2:32];
Rows=[StartRow:2:16];

PatternsToUse=[];
for i=1:length(Patterns)
    Electrode=Patterns(i);
    Column=Indexes(Electrode,1);
    Row=Indexes(Electrode,2);
    a1=find(Columns==Column)
    a2=find(Rows==Row)
    if (a1 & a2)
        PatternsToUse=[PatternsToUse Electrode];
    end
end

Electrodes=NS512_PatternsForRectangularAreaStimulation(RowsIndexes,ColumnsIndexes);

Chunks=[];
for i=1:8
    PatternsForMovie=Patterns((i-1)*16+1:i*16);
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end

MovieChunksFile=[8 Chunks];

break;

cd F:\StimFiles; 
fid = fopen('128el_30ms_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('128el_30ms_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('128el_30ms_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
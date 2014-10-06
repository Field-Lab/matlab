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

Columns1=[1:4:32];
Rows1=[1:4:16];

Columns2=[3:4:32];
Rows2=[3:4:16];

PatternsToUse1=[];
PatternsToUse2=[];
for i=1:length(Patterns)
    Electrode=Patterns(i);
    Column=Indexes(Electrode,1);
    Row=Indexes(Electrode,2);
    a1=find(Columns1==Column)
    a2=find(Rows1==Row)
    if (a1 & a2) 
        PatternsToUse1=[PatternsToUse1 Electrode];
    end
    a1=find(Columns2==Column)
    a2=find(Rows2==Row)
    if (a1 & a2)
        PatternsToUse2=[PatternsToUse2 Electrode];
    end
end

Chunks=[];
for i=1:4
    PatternForMovie1=PatternsToUse1((i-1)*8+1:i*8);
    P1=[PatternForMovie1 PatternForMovie1([4:-1:1 8:-1:5])]
    PatternForMovie2=PatternsToUse2((i-1)*8+1:i*8);
    P2=[PatternForMovie2 PatternForMovie2([4:-1:1 8:-1:5])]
    %PatternsForMovie=Patterns((i-1)*16+1:i*16);
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    Chunk1=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,P1);
    Chunk2=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,P2);
    Chunks=[Chunks Chunk1 Chunk2];
end

MovieChunksFile=[8 Chunks];
%break     
cd F:\StimFiles; 
fid = fopen('64el_30ms_2seq_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('64el_30ms_2seq_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('64el_30ms_2seq_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
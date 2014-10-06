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
Columns=Indexes(:,1);
Rows=Indexes(:,2);

StartColumn=[8 6 8];
ColumnsNumber=[5 6 5];
ElectrodesToUse=[];

for i=1:3
    RowID=4+(i-1)*4;
    for j=1:ColumnsNumber(i)
        ColumnID=StartColumn(i)+(j-1)*4;
           
        for Pattern=Patterns
            a1=find(Columns(Pattern)==ColumnID);
            a2=find(Rows(Pattern)==RowID);            
            if (a1 & a2)
                ElectrodesToUse=[ElectrodesToUse Pattern];
            end
        end
    end
end

NumberOfChunks=120;
Amplitudes=[0.40 0.74 1.36 2.50];
Array=zeros(512,1920);
PatternID=0;
for i=1:16
    El1=ElectrodesToUse(i);
    for j=i+1:16        
        El2=ElectrodesToUse(j);
        for amp1=1:4
            for amp2=1:4
                PatternID=PatternID+1;
                Array(El1,PatternID)=Amplitudes(amp1);
                Array(El2,PatternID)=Amplitudes(amp2);
            end
        end
    end
end

Chunks=[];
for i=1:NumberOfChunks
    PatternsForMovie=[1:16]+(i-1)*16;
    %PatternsForMovie=Patterns((i-1)*16+1:i*16);
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end

MovieChunksFile=[NumberOfchunks Chunks];

cd C:\home\Pawel\nauka\StimFiles\2012_02_02; 
fid = fopen('16el_all_pairs_30ms_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('16el_all_pairs_30ms_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('16el_all_pairs_30ms_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
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

StartColumn=[4 5 6 5 4];
ColumnsNumber=[13 13 12 13 13];
PatternsToUse=[];

for i=1:5
    RowID=4+(i-1)*2;
    for j=1:ColumnsNumber(i)
        ColumnID=StartColumn(i)+(j-1)*2;
           
        for Pattern=Patterns
            a1=find(Columns(Pattern)==ColumnID);
            a2=find(Rows(Pattern)==RowID);            
            if (a1 & a2)
                PatternsToUse=[PatternsToUse Pattern];
            end
        end
    end
end

Chunks=[];
for i=1:4
    PatternsForMovie=PatternsToUse((i-1)*16+1:i*16);
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end

MovieChunksFile=[4 Chunks];
break
cd C:\home\Pawel\nauka\StimFiles\2012_02_02; 
fid = fopen('64el_5_rows_30ms_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('64el_5_rows_30ms_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('64el_5_rows_30ms_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
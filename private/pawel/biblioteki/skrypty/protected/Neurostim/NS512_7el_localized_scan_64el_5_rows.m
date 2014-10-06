electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

TimeShiftInMs=30;
DelayInMs=30;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=zeros(512,64); 
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

for i=1:length(PatternsToUse)
    Radius=1;
    AllElectrodes=electrodeMap.getAdjacentsTo(PatternsToUse(i),Radius)';
    
    Array(PatternsToUse(i),i)=1;
    for j=2:length(AllElectrodes)
        Array(AllElectrodes(j),i)=0.16667;
    end
end


Chunks=[];
for i=1:4
    %PatternsForMovie=Patterns((i-1)*16+1:i*16);
    PatternsForMovie=[1:16]+(i-1)*16;
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end
break
MovieChunksFile=[4 Chunks];

cd C:\home\Pawel\nauka\StimFiles\2012_02_02; 
fid = fopen('64el_7el_per_pattern_5_rows_30ms_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('64el_7el_per_pattern_5_rows_30ms_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('64el_7el_per_pattern_5_rows_30ms_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
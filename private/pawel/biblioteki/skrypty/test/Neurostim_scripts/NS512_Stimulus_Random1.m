TimeShiftInMs=30;
DelayInMs=30;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
%Array=zeros(512,64); 
%Array=eye(512,512); 
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


NumberOfMovies=20;
NumberOfElectrodes=length(ElectrodesToUse);
ElectrodesCoeffs=zeros(1,NumberOfElectrodes);
Array=[];
Counter=0;
Chunks=[];
for i=1:NumberOfMovies
    PatternsForMovie=[];
    TimesForMovie=[];
    for Pattern=1:16        
        TimesOffset=TimeShift+Delay*(Pattern-1);
        
        %Electrodes=ceil(rand(1,4)*NumberOfElectrodes);
        Electrodes=NS512_UniqueRandomIntNumbers(NumberOfElectrodes,4);
        Times=ceil(rand(1,4)*10);
    
        UniqueTimes=unique(Times);
        for TimeID=UniqueTimes
            Counter=Counter+1;
            ElectrodesCoeffs=zeros(512,1);
            WhichElectrodes=find(Times==TimeID);
            ElectrodesCoeffs(ElectrodesToUse(Electrodes(WhichElectrodes)))=1;
            Array=[Array ElectrodesCoeffs];
            PatternsForMovie=[PatternsForMovie Counter];
            TimesForMovie=[TimesForMovie TimesOffset+TimeID];
            %Times(=TimesOffset+TimeID;
        end
    end 
    Chunk=NS_MovieChunkGenerationForExperiment(TimesForMovie,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end

MovieChunksFile=[NumberOfMovies Chunks];

cd C:\home\Pawel\nauka\StimFiles\2012_02_02; 
fid = fopen('16el_random1_30ms_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('16el_random1_30ms_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('16el_random1_30ms_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
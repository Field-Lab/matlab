%javaaddpath('C:\pawel\nauka\Vision\vision.jar')

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

Pairs=NS512_PairsFor16Electrodes;

%test:
%%figure(1);
%clf
%hold on;
%c=zeros(16,16);
%for i=1:15
%    for j=1:8
%        x=Pairs(i,j,1);
%        y=Pairs(i,j,2);
%        c(x,y)=c(x,y)+1;
%        if x>y
%            plot(x,y,'bd');
%        else
%           plot(y,x,'bd');
%       end
%        axis([0 17 0 17]);
%    end
%nd
   
AllAmpCombinations=zeros(16,2);
CombID=0;
for i=1:4
    for j=1:4
        CombID=CombID+1;
        AllAmpCombinations(CombID,1)=Amplitudes(i);
        AllAmpCombinations(CombID,2)=Amplitudes(j);
    end
end

Chunks=[];
Times=[TimeShift:Delay:TimeShift+15*Delay];
for P=1:8 %8 chunks will include stimulation for a specific pair of electrodes (in given chunk always two combinations of amplitudes)
    AmpCombinations=[(P-1)*2+1 (P-1)*2+2]; 
    for Round=1:15 %one round is one movie chunk
        ChunkID=(P-1)*15+Round
        for PairID=1:8            
            PatternID=(ChunkID-1)*16+PairID
            El1=Pairs(Round,PairID,1);
            El2=Pairs(Round,PairID,2);
            Array(El1,PatternID)=AllAmpCombinations(AmpCombinations(1),1);
            Array(El2,PatternID)=AllAmpCombinations(AmpCombinations(1),2);
            Array(El1,PatternID+8)=AllAmpCombinations(AmpCombinations(2),1);
            Array(El2,PatternID+8)=AllAmpCombinations(AmpCombinations(2),2);
            PatternsForMovie(PairID)=PatternID;
            PatternsForMovie(PairID+8)=PatternID+8;            
        end
        Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
        Chunks=[Chunks Chunk];
    end
end

%test2        
%for i=1:16
%    El1=ElectrodesToUse(i);
%    for j=i+1:16        
%        El2=ElectrodesToUse(j);
%        for amp1=1:4
%            for amp2=1:4                
%                a=length(find(Array(El1,:)==Amplitudes(amp1) & Array(El2,:)==Amplitudes(amp2)));
%                if length(a)~=1
%                    'dfgdfgdfgdfgdfg'
%                    El1
%                    El2
%                    Amplitudes(amp1)
%                    Amplitudes(amp2)
%                    a
%                end                
%            end
%        end
%    end
%end
%break
MovieChunksFile=[ChunkID Chunks];

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
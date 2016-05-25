TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

Patterns=NS512_OptimalElectrodeSequence();

Times=[TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];

cd /Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_files_test; 
           
Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);
Columns=Indexes(:,1);

AllPatterns=[];

for i=1:4
    Chunks=[];
    PatternsForSet=[];
    Column_start=(i-1)*8;
    Cols=find(Columns>Column_start & Columns<Column_start+9); %128 electrodes to be used now
    for j=1:512
        if find(Cols==Patterns(j))
            PatternsForSet=[PatternsForSet Patterns(j)]; % the optimal sequece of these 128 electrodes
        end            
    end
        
    for j=1:2
        Start=(j-1)*64;
        PatternsForMovie=PatternsForSet(Start+1:Start+64);
        Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
        Chunks=[Chunks Chunk];
        AllPatterns=[AllPatterns PatternsForMovie]; % just to check
    end
    MovieChunksFile=[2 Chunks];
    name=['128el_mv' num2str(i)];
    fid=fopen(name,'wb');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid); 
end

fid = fopen('128el_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('128el_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

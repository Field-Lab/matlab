%This scripts generates two sequenecs of electrodes: 
%1) 128 electrodes - every second row and pattern
%20 128 electrodes - all the electrodes in the center quarter of the array
%(8x16)
%Then, the script generates three sets of stimulus files: with the first
%sequence, the second one, and with combined sequence (224 electrodes).

TimeShiftInMs=30;
DelayInMs=30;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

StartRow=1;
StartColumn=1;
ElectrodesSparse=NS512_PatternsForRectangularAreaStimulation([StartRow:2:16],[StartColumn:2:32]);
ElectrodesCenter=NS512_PatternsForRectangularAreaStimulation([5:12],[9:24]);

ElectrodeOrder=NS512_OptimalElectrodeSequence();
ElectrodesSparseOrder=[];
ElectrodesCenterOrder=[];
ElectrodesCombinedOrder=[];

for i=ElectrodeOrder
    if find(ElectrodesSparse==i)
        ElectrodesSparseOrder=[ElectrodesSparseOrder i];
    end
    if find(ElectrodesCenter==i)
        ElectrodesCenterOrder=[ElectrodesCenterOrder i];
    end
    if (find(ElectrodesSparse==i |ElectrodesCenter==i))
        ElectrodesCombinedOrder=[ElectrodesCombinedOrder i];
    end
end

ChunksSparse=[];
ChunksCenter=[];
for i=1:8
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    
    ElectrodesSparseForMovie=ElectrodesSparseOrder((i-1)*16+1:i*16);
    ChunkSparse=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,ElectrodesSparseForMovie);
    ChunksSparse=[ChunksSparse ChunkSparse];
    
    ElectrodesCenterForMovie=ElectrodesCenterOrder((i-1)*16+1:i*16);
    ChunkCenter=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,ElectrodesCenterForMovie);
    ChunksCenter=[ChunksCenter ChunkCenter];
end

ChunksCombined=[];
for i=1:14
    Times=[TimeShift:Delay:TimeShift+15*Delay];
    
    ElectrodesCombinedForMovie=ElectrodesCombinedOrder((i-1)*16+1:i*16);
    ChunkCombined=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,ElectrodesCombinedForMovie);
    ChunksCombined=[ChunksCombined ChunkCombined];
end

MovieChunksFileSparse=[8 ChunksSparse];
MovieChunksFileCenter=[8 ChunksCenter];
MovieChunksFileCombined=[14 ChunksCombined];

cd C:\home\pawel\nauka\brain_slices\StimulusFiles;
fid = fopen('All_el','wb');
fwrite(fid,electrodes,'int32');
fclose(fid);
fid = fopen('All_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);
fid = fopen('Sparse_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFileSparse,'int32');
fclose(fid); 
fid = fopen('Center_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFileCenter,'int32');
fclose(fid); 
fid = fopen('Combined_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFileCombined,'int32');
fclose(fid); 
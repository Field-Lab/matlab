function MovieChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);

Array=zeros(1,6+length(Patterns)*3);
Array(1,6)=NumberOfSamples;

for i=1:length(Patterns)
    index=6+(i-1)*3;
    Array(1,index+1)=Times(i); %TimeShift+(i-1)*Delay; %point in time for this pattern
    Array(1,index+2)=Patterns(i);
    Array(1,index+3)=1; %scaling factor - for now always equal to 1
end

MovieChunk=[length(Array) Array];
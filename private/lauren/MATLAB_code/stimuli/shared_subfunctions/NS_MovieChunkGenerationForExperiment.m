function MovieChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns)

%arguments
%   Times: vector of time values (in samples) at which corresponding patterns will be played
%   NumberOfSamples: how long the total movie chunk lasts in samples (check)
%   Patterns: vector of pattern numbers that corresponds to particular patterns stored in patterns
%   stimulus files (must be same length as Times)

Array=zeros(1,6+length(Patterns)*3); %first 5 values after Array length see below) will be left at 0
Array(1,6)=NumberOfSamples; %7th value specifies length of movie chunk in samples

for i=1:length(Patterns)
    index=6+(i-1)*3;
    Array(1,index+1)=Times(i); %TimeShift+(i-1)*Delay; %point in time for this pattern
    Array(1,index+2)=Patterns(i); %pattern to be played
    Array(1,index+3)=1; %scaling factor - for now always equal to 1
end

MovieChunk=[length(Array) Array]; %first value specifies number of elements in rest of movieChunk
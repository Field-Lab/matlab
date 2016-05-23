TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;

TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 

Patterns=[1:519];
Disconnected=[1 130 259 260 389 390 519];
PatternsToUse=NS_RemoveBadChannels(Patterns,Disconnected); %at this point we have the electrode order we want; this
%electrode sequence includes all the electrodes that are actually connected
%to the electronics (there are of course 512 of them).
% Please redefine this sequence up to your taste.

%Now apply the electrodes shuffle
PatternsMessed=NS512_519InverseMess(PatternsToUse);

% And now just generate the stim files
Times=[TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];

Chunks=[];
p=[];
for i=1:8
    PatternsForMovie=PatternsMessed((i-1)*64+1:i*64)
    %p=[p PatternsForMovie]
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end
%break;
MovieChunksFile=[i Chunks];
%break;
keyboard; 
% cd C:\pawel\nauka\granty\NCN_Harmonia2013\realizacja\stim_files\scan519; 
fid = fopen('519_512_el','wb')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('519_512_pt','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('519_512_mv','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
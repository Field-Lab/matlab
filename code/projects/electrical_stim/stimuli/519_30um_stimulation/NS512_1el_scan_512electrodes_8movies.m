% Script to generate a single electrode scan on the 30um stimulation board.

% Define stimulation time parameters
TimeShiftInMs=7.5;
InterPulseLatencyInMs=7.5;
TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;
NumberOfSamples=10000;

% Hack - the software still assumes that the 512-stim board is the only one
% that exists for electrical stimulation. Therefore, the electrodes file
% (and the Array matrix) can only contain 512 electrodes. 
Electrodes = 1:512;
Array = eye(512,512); 

% Pattern numbers correspond to the 519 electrode map. 
Patterns = 1:519;

% Since 519>512, some channels are disconnected
Disconnected=[1 130 259 260 389 390 519];
PatternsToUse=NS_RemoveBadChannels(Patterns,Disconnected);

% at this point we have the electrode order we want; this
%electrode sequence includes all the electrodes that are actually connected
%to the electronics (there are of course 512 of them).
% Please redefine this sequence up to your taste.

% Now apply the electrode transform to map the desired 519 patterns on the 
% 512 system
PatternsTransformed=NS512_519Inverse(PatternsToUse);

% And now just generate the stim files
Times = TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency;

Chunks=[];
p=[];
for ii=1:8
    PatternsForMovie = PatternsTransformed((ii-1)*64+1:ii*64);
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    Chunks=[Chunks Chunk];
end
MovieChunksFile=[ii Chunks];

keyboard; 

fid = fopen('519_512_el','wb');
fwrite(fid,Electrodes,'int32');
fclose(fid);

fid = fopen('519_512_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('519_512_mv','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
%% Script to generate the full array single electrode scan with local return

% Set the size of the trigger interval in number of samples
numberOfSamples = 10000; %trigger interval = 0.5 s (standard for STIM512)
% Set an offset time from the beginning of the stimulus
timeShiftInMs = 7.5;
% Set the desired time between pulses
interPulseLatencyInMs = 7.5;

% Convert to sample points
timeShift=timeShiftInMs*20;
interPulseLatency=interPulseLatencyInMs*20;

% Create electrodes vector
electrodes = 1:512;

% Create array variable to define stimulus patterns
array = eye(512,512); 
%Add local returns
for e = 1:512
    cluster = getCluster512(e); 
    array(e,cluster(2:end)) = -1/(length(cluster)-1); 
end

% Get sequence of stimulus patterns
patterns = NS512_OptimalElectrodeSequence();

% Create time vector to present each pulse 63 hard coded?
t = timeShift:interPulseLatency:timeShift+63*interPulseLatency;
% cd /Users/grosberg/Desktop/  
%%
allPatterns=[];
for i=1:8 % 8 Movie ChunksA
    Start=(i-1)*64;
    patternsForMovie=patterns(Start+1:Start+64);
    mChunk=NS_MovieChunkGenerationForExperiment(t,numberOfSamples,patternsForMovie);
    allPatterns=[allPatterns mChunk];
end
% 8 movie chunks. 
MovieChunksFile=[8 allPatterns];
%%
keyboard;
%%
fid=fopen('512el_lr_mv','wb');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);

fid = fopen('512el_lr_el','wb');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('512el_lr_pt','wb','ieee-le.l64');
fwrite(fid,array,'double');
fclose(fid);

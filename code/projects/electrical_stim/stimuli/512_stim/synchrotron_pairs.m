% two random orders in 2 movie chunks?

% Stimulates at the center electrode and all nearest neighbor pairs (6) as
% well as all pairs of nearest neighbors (6)

% Delays the start of the first stimulus.  Almost always set to 0.
TimeShiftInMs=0;

% Center electrodes are picked to stimulate axons.  Pick low threshold
% electrodes ideally with a neighboring electrode also stimulating the
% axon. Good distance from initial segment.

CenterElectrodes=[367 90]; %limited to 2 centers % 367 targets ON parasol 421, and 1265
% CenterElectrodes=[363 282]; %limited to 2 centers % 367 targets ON parasol 421, and OFF parasol 4218
% CenterElectrodes=[21 93]; %limited to 2 centers % 21 targets ON parasol 124, and 93 targets OFF parasol 1381

NumberOfClusters=length(CenterElectrodes);

DelayInMs = 30; %set trigger to 2 s

timeShiftInMs      = 0;        % offset the first stimulus from the beginning movie chunk
arraySpacing       = 60;       % 30 or 60 micron spacing
delayInMs          = 40;       %    
interPairDelayInMs = 7.5;       % clusterDelayInMs*numClusters < delayInMs
numberOfSamples    = 10000;    % 0.5s SET TO TRIGGER INTERVAL IN LABVIEW

%%
sameAmps = 1; 
[electrodes Array] = generatePatternsSynchrotronStim(CenterElectrodes, sameAmps);

%% Generate movie files
%convert inputs to samples
Delay          = delayInMs*20;
TimeShift      = timeShiftInMs*20; % in sampling periods (50 microseconds)
interPairDelay = interPairDelayInMs*20; 
pulseDelay1    = 10; % 0.5 ms delay
pulseDelay2    = 20; % 100us delay
pulseDelay3    = 40; % 1 ms delay

% Times for zero time spacing
Times = (0:interPairDelay:interPairDelay*(48-1)); 
Times1_1 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48;
Times1 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48 + pulseDelay1; 
Times2_1 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48 + interPairDelay*24 ;
Times2 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48 + interPairDelay*24 + pulseDelay2; 
Times3_1 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48 + interPairDelay*24*2 ; 
Times3 = (0:interPairDelay:interPairDelay*(24-1)) + interPairDelay*48 + interPairDelay*24*2  + pulseDelay3; 

allTimes = sort(cat(2,Times,Times1,Times1_1,Times2,Times2_1,Times3_1,Times3)); 


Times1 = allTimes(1:76); 
Patterns1 = [1:1:24 ;39:1:62]; 
firstPatterns= Patterns1(:); 
%% horrible hand coding
i=1; 
for e = 1:6
    order(i) = 25; 
    i = i + 1; 
    secondE = 25 + 7 + e; 
    order(i) = secondE; 
    i = i+1; 
end
for e = 1:6
    order(i) = 25 + e; 
    i = i + 1; 
    order(i) = 32;
    i = i+1; 
end

% group = [32 34:38];
% for e = 1:6
%     order(i) = 26; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
% group = [32:33 35:38];
% for e = 1:6
%     order(i) = 27; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
% group = [32:34 36:38];
% for e = 1:6
%     order(i) = 28; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
% group = [32:35 37:38];
% for e = 1:6
%     order(i) = 29; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
% group = [32:36 38];
% for e = 1:6
%     order(i) = 30; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
% group = [32:37];
% for e = 1:6
%     order(i) = 31; 
%     i = i + 1; 
%     order(i) = group(e); 
%     i = i+1; 
% end
%%
order2 = order + 38; 
%%
order = reshape(order,2,[]); 
order2 = reshape(order2,2,[]); 
a = order';
b = order2';
col_interleave = reshape([a(:) b(:)]',2*size(a,1), [])'; 
allOrders = col_interleave(:);

%%
patterns = cat(1, firstPatterns,repmat(allOrders,3,1)); 
%% how many chunks?
nChunksTotal = ceil(allTimes(end)/10000); 
movieChunks = nChunksTotal;

startIndex = 1; 
figure; 
for i = 1:nChunksTotal
    lastTimePoint = 10000*i;
    endIndex = find(allTimes<lastTimePoint,1,'last'); 
    thisPatternChunk = patterns(startIndex:endIndex)'; 
    Times = allTimes(startIndex:endIndex)-10000*(i-1);
    hold on; plot(Times,1,'-x');
    startIndex = endIndex+1; 
    Chunk=NS_MovieChunkGenerationForExperiment(Times, numberOfSamples, thisPatternChunk);
    movieChunks = [movieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

MovieChunksFile=movieChunks;

%%
keyboard;

fid = fopen('synch_sameAmps_el','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('synch_sameAmps_pt','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('synch_sameAmps_mv','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
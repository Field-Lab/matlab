% Script to generate pulse triplets


%% Test 1: single electrode

TimeShiftInMs = 0;
DelayInMs     = 7.5; %time between subsequent pulses on different electrodes

NumberOfSamples = 10000;     % 2s SET TO TRIGGER INTERVAL IN LABVIEW
electrodes = [1 3 6 7]; % Choose 4 electrodes
% NxM N - number of electrodes used in all the movie files
%     M - number of different patterns in all the movies in a given run

C = 2.5;  % Maximum stimulus current applied to any one electrode in the pattern
nPatterns = 400000; 
Array = zeros(size(electrodes,2),nPatterns); 
for p = 1:nPatterns
Array(:,p) = rand(1,size(electrodes,2))*2*C - C; % Randomly assign stimulus amplitudes on [-C,C] to electrodes
end

TimeShift   = TimeShiftInMs*20;
Delay       = DelayInMs*20;
Times       = (0:Delay:(NumberOfSamples-Delay-TimeShift)) + TimeShift; 
MovieChunks = ceil(nPatterns/size(Times,2));

disp(['record for ' num2str(MovieChunks*NumberOfSamples/20000) ' s']); 

length      = size(Times,2);
Pattern     = [1:1:nPatterns 1:(length*MovieChunks-nPatterns)]; 
Pattern     = reshape(Pattern, MovieChunks, []);

for i =1:MovieChunks
    Patterns = Pattern(i,1:length);
    Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end
%%
keyboard; 

% Write files
cd ~/Desktop/StimFiles
fid = fopen([num2str(nPatterns) 'Patterns_electrodes'],'wb')
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen([num2str(nPatterns) 'Patterns_patterns'],'wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen([num2str(nPatterns) 'Patterns_movie'],'wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 


TimeShiftInMs           = 7.5;
InterPulse1LatencyInMs  = 5;
InterPulse2LatencyInMs  = 1.5; %1.5?
InterTripletLatencyInMs = 7.5;

TimeShift               = TimeShiftInMs*20;
InterPulse1Latency      = InterPulse1LatencyInMs*20;
InterPulse2Latency      = InterPulse2LatencyInMs*20; %1.5?
InterTripletLatency     = InterTripletLatencyInMs*20;
% InterPulseLatency       = InterPulseLatencyInMs*20;

NumberOfSamples         = 10000; %trigger interval = 0.5 s
electrodes=[1:512];
Array=eye(512,512); 

Array    = repmat(Array,1,7);

Patterns = NS512_OptimalElectrodeSequence();

% Times   = [TimeShift:InterPulseLatency:TimeShift+63*InterPulseLatency];
repPeriod = InterPulse1Latency + InterPulse2Latency + InterTripletLatency; 

numSequences = floor((NumberOfSamples-TimeShift)/repPeriod);

pulse1_times = TimeShift + 0:repPeriod:numSequences*repPeriod; 
pulse2_times = TimeShift + InterPulse1Latency + 0:repPeriod:numSequences*repPeriod; 
pulse3_times = TimeShift + InterPulse1Latency + InterPulse2Latency + 0:repPeriod:numSequences*repPeriod; 

Times = sort(cat(2,pulse1_times,pulse2_times,pulse3_times)); 
% cd /Volumes/Stream-phoenix/Analysis/stim512/2012-09-18-1/stim_files_test; 
cd /Users/grosberg/Desktop/  

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

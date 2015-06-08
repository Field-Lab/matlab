% Script to generate pulse triplets for artifact recordings
% L Grosberg 7/2014

%% Test 1: single electrode stimulation

TimeShiftInMs           = 0;
InterPulse1LatencyInMs  = 20;  % time between first two pulses
InterPulse2LatencyInMs  = 1.5; % time between second and third pulses
TripletRepInMs          = 30;

TimeShift               = TimeShiftInMs*20;
InterPulse1Latency      = InterPulse1LatencyInMs*20;
InterPulse2Latency      = InterPulse2LatencyInMs*20; %1.5?
TripletRep              = TripletRepInMs*20;

NumberOfSamples         = 10000; % trigger interval = 0.5 s
electrodes              = 1:512;
Array                   = eye(512,512);
prePulseAmp             = 4; % prepulse at 4uA

% ampsToTest              = [0.189  0.2835  0.378   0.4725  0.567];
% ampsToTest              = [0.6615 0.756   0.8505  0.945   1.0395];
% ampsToTest              = [1.134  1.2285  1.323   1.512   1.6065];
 ampsToTest              = [1.7955 2       2.1735  2.3625  2.646];

numAmpsPerElec          = size(ampsToTest,2);    
Array                   = repmat(Array,1,2+numAmpsPerElec); % 2 patterns/electrode associated with first two large pulses, numAmpsPerElec associated with the amps to be tested
% NxM N - number of electrodes used in all the movie files
%     M - number of different patterns in all the movies in a given run

%%
multiplier   = [prePulseAmp prePulseAmp ampsToTest]; % Two large amplitude pre-pulses and the amplitudes to test
multiplier   = repmat(reshape(repmat(multiplier,512,1),1,[]),512,1);
Array        = Array.*multiplier; 
numSequences = floor((NumberOfSamples-TimeShift)/TripletRep);

pulse1_times = TimeShift + 0:TripletRep:numSequences*(TripletRep-1); 
pulse2_times = TimeShift + InterPulse1Latency + 0:TripletRep:numSequences*(TripletRep-1); 
pulse3_times = TimeShift + InterPulse1Latency + InterPulse2Latency + 0:TripletRep:numSequences*(TripletRep-1); 

Times        = sort(cat(2,pulse1_times,pulse2_times,pulse3_times)); 

Patterns     = NS512_OptimalElectrodeSequence();
Patterns     = repmat(Patterns,[3*numAmpsPerElec 1]); 

% addM         = (0:(size(Patterns,1)-1))*512; 
% addM         = repmat(addM',[1 512]);
% Patterns     = Patterns + addM; 
%%
test     = [zeros(numAmpsPerElec,1) ones(numAmpsPerElec,1) (2:(numAmpsPerElec+1))']';
test     = test(:).*512; 
testM    = repmat(test,[1 512]); 
Patterns = Patterns + testM; 
 
pat = [];
for i    = 1:numAmpsPerElec
    temp = Patterns(3*(i-1)+(1:3),:);
    pat  = cat(1,pat,temp(:));
end

pat          = pat'; 
nPatterns    = size(pat,2); 
len          = size(Times,2);
nChunksTotal = ceil(nPatterns/len);

    
Patterns     = [pat pat(1:(len*nChunksTotal-nPatterns))];
Patterns     = reshape(Patterns, [], nChunksTotal);
%%
movieChunks  = nChunksTotal; 
for i =1:nChunksTotal
    currentPattern = Patterns(1:len,i)';
    Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, currentPattern);
    movieChunks = [movieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end
% %%
%  cd /Users/grosberg/Desktop/  
% 
% fid = fopen('pulse_triplets_4uAampSet4_el','wb','ieee-le.l64');
% fwrite(fid,electrodes,'int32');
% fclose(fid);
% 
% fid = fopen('pulse_triplets_4uAampSet4_pt','wb','ieee-le.l64');
% fwrite(fid,Array,'double');  
% fclose(fid);     
% 
% fid = fopen('pulse_triplets_4uAampSet4_mv','wb','ieee-le.l64');
% fwrite(fid,movieChunks,'int32');
% fclose(fid);
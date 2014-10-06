function MovieChunksFile = generateMovieClusterStim2SecondariesTest(nClusters,TimeShiftInMs,DelayInMs,nSamples)

% generates movieChunksFile for stimulus file generation
%
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% movie chunk length is limited to <=2 seconds, so patterns must be broken up into 11 movie chunks, corresponding
% to different rows in Pattern
%
% *** patterns must be for full 7-electrode clusters (265 patterns/cluster) because pattern number
% are hard-coded in this function! 
%
% order of patterns played:
% movie chunk 1:
%   first half (1:nClusters*24)
%   random primary + pair of secondaries (1:0.5:0.5) from cluster 1 (chosen from all polarity combinations and all
%   secondary electrode pairs; 24 total to choose from)
%   random primary + pair of secondaries from cluster 2
%   random primary + pair of secondaries from cluster 3
%   different random primary + pair of secondaries from cluster 1 (then cluster 2, then cluster 3)
%   etc.
%
%   second half (nCluster*24+1:nClusters*48) same as first, but new random orders
%
% movie chunk 2: as chunk 1, but with 1:1:0.5 amplitude ratio
% movie chunk 3: as chunk 1, but with 1:2:0.5 amplitude ratio
% movie chunk 4: as chunk 1, but with 1:0.5:1 amplitude ratio
% movie chunk 5: as chunk 1, but with 1:1:1 amplitude ratio
% movie chunk 6: as chunk 1, but with 1:2:1 amplitude ratio
% movie chunk 7: as chunk 1, but with 1:0.5:2 amplitude ratio
% movie chunk 8: as chunk 1, but with 1:1:2 amplitude ratio
% movie chunk 9: as chunk 1, but with 1:2:2 amplitude ratio

% movie chunk 10:
%   first quarter (1:nClusters*12): 
%   random primary + secondary from cluster 1 with 1:0.5 amplitude ratio (both polarities on secondary)
%   random primary + secondary from cluster 2 with 1:0.5 amplitude ratio
%   random primary + secondary from cluster 3 with 1:0.5 amplitude ratio
%
%   second quarter: as first, but new random order
%   third quarter: as first, but 1:1 amplitude ratio
%   fourth quarter: as third, but new random order
%   
% movie chunk 11: 
%   first quarter: random primary + secondary from cluster 1 with 1:2 amplitude ratio
%   second quarter: as first quarter, but new random order
%   third quarter: random individual electrode (primary positive or secondary either polarity)
%   fourth quarter: as third, but new random order
%
% author: Lauren (SNL-E), based on function written by Pawel Hottowy
% 2009-08-24

Pattern = zeros(1, 50*nClusters);
PatternsPerCluster = 90; %number of patterns in one cluster

keyboard
% values in Pattern identify order in wich stimulus patterns (specified in patterns array, which is generated
% by a different function) are played
for i=1:nClusters
    offset = PatternsPerCluster*(i-1);
    
    % primary plus pair of neighboring electrodes
    Pattern(1,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 36; % center with 0.5, 0.5 secondaries
    Pattern(2,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 60; % center with 1, 0.5 secondaries
    Pattern(3,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 84; % center with 2, 0.5 secondaries
    Pattern(4,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 108; % center with 0.5, 1 secondaries
    Pattern(5,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 132; % center with 1, 1  secondaries   
    Pattern(6,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 156; % center with 2, 1 secondaries
    Pattern(7,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 180; % center with 0.5 2 secondaries
    Pattern(8,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 204; % center with 1 2 secondaries
    Pattern(9,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 228; % center with 2 2 secondaries
    
    %primary plus individual secondaries
    Pattern(10,i:nClusters:i+nClusters*47) = [randperm(12) randperm(12) randperm(12)+12 randperm(12)+12] + offset; % center with single secondary (0.5x and 1x)
    
    %primary plus individual secondaries and each electrode individually
    Pattern(11,i:nClusters:i+nClusters*49) = [randperm(12)+24 randperm(12)+24 randperm(13)+252 randperm(13)+252] + offset; % center with single secondary (2x) and each electrode individually
end

% checks to make sure that each pattern number is in Pattern twice--for initial testing
% for i = 1:795
%     if sum(sum(Pattern == i)) ~= 2
%         disp(['pattern' num2str(i) 'is not in Pattern twice'])
%     end
% end

TimeShift = TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay = round(DelayInMs*20/nClusters); %number of samples between patterns on different clusters

% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = 11;
for i=1:11
    if i<11
        length = 48*nClusters;
    else %movie chunk that includes primary-alone patterns
        length = 50*nClusters;
    end
    Patterns = Pattern(i, 1:length);
    Times = TimeShift : Delay : TimeShift+Delay*(length-1);
    Chunk = NS_MovieChunkGenerationForExperiment(Times, nSamples, Patterns);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

MovieChunksFile=MovieChunks;


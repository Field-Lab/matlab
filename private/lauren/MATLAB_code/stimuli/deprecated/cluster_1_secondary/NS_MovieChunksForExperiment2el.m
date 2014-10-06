function [PatternOrder, MovieChunksFile]=NS_MovieChunksForExperiment2el(nClusters,TimeShiftInMs,DelayInMs,nSamples)

% generates movieChunksFile for stimulus file generation
%
% specifies order in which patterns are played during stimulus (time each pattern is played)
%
% movie chunk length is limited so patterns must be broken up into 5 movie chunks, corresponding
% to different rows in Pattern (check)
%
% *** patterns must be for full 7-electrode clusters (110 patterns/cluster) because pattern number
% are hard-coded in this function! 
%
% order of patterns played:
% movie chunk 1:
%   first half (1:nClusters*24)
%   random primary + 0.25 secondary from cluster 1 (chosen from all polarity combinations and all
%   secondary electrodes; 24 total to choose from)
%   random primary + 0.25 secondary from cluster 2
%   random primary + 0.25 secondary from cluster 3
%   different random primary pos. + 0.25 secondary pos. from cluster 1 (then cluster 2, then cluster 3)
%   etc.
%
%   second half (nCluster*24+1:nClusters*48) same as first, but new random orders
%
% movie chunk 2: same as movie chunk 1, but amplitude ratio is 1 primary : 0.5 secondary
%
% movie chunk 3: first pattern in each half corresponds to center electrode alone, and the rest of 
% the patterns in each half are the same as movie chunk 1, but amplitude ratio is 1 primary : 0.75 secondary
%
% movie chunk 4: same as chunk 3, but with negative polarity for center electrode alone
%
% movie chunk 5: 
%   each quarter:
%   random secondary electrode alone from cluster 1 (either polarity)
%   random secondary electrode alone from cluster 2 (either polarity)
%   random secondary electrode alone from cluster 3 (either polarity)
%   different random secondary electrode alone from cluster 1 (either polarity)
%   etc.


Pattern = zeros(5,50*nClusters);
PatternsPerCluster = 110; %number of patterns in one cluster, with all polarity combinations,
% 4 primary:secondary amplitude ratios, and each electrode stimulated individually (both polarities)


% values in Pattern identify order in wich stimulus patterns (specified in patterns array, which is generated
% by a different function) are played
for i=1:nClusters
    offset = PatternsPerCluster*(i-1);
    Pattern(1,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset; % center with 0.25 for secondary el.
    Pattern(2,i:nClusters:i+nClusters*47) = [randperm(24) randperm(24)] + offset + 24; % center with 0.5 for secondary el.
    Pattern(3,i:nClusters:i+nClusters*49) = [97+offset randperm(24)+offset+48 97+offset randperm(24)+offset+48]; % center with 0.75 for secondary el., plus center alone positive
    Pattern(4,i:nClusters:i+nClusters*49) = [98+offset randperm(24)+offset+72 98+offset randperm(24)+offset+72]; % center with 1 for secondary el., plus center alone negative
    Pattern(5,i:nClusters:i+nClusters*47) = [randperm(12) randperm(12) randperm(12) randperm(12)] + offset + 98; % only secondary electrodes (both polarities)
end

% from original version of function (only 98 patterns/cluster)
% Pattern=zeros(nClusters,50*nClusters);
% for i=1:nClusters
%     Pattern(1,i:nClusters:i+nClusters*47)=[randperm(24) randperm(24)]+98*(i-1); 
%     Pattern(2,i:nClusters:i+nClusters*47)=[randperm(24) randperm(24)]+98*(i-1)+24;
%     Pattern(3,i:nClusters:i+nClusters*49)=[97+98*(i-1) randperm(24)+98*(i-1)+48 97+98*(i-1) randperm(24)+98*(i-1)+48];
%     Pattern(4,i:nClusters:i+nClusters*49)=[98+98*(i-1) randperm(24)+98*(i-1)+72 98+98*(i-1) randperm(24)+98*(i-1)+72];
% end

TimeShift = TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay = round(DelayInMs*20/nClusters); %number of samples between patterns on different clusters

% makes a movie chunk for each row of Pattern and concatenates MovieChunks
MovieChunks = 5;
for i=1:5
    if i<3 || i==5
        length = 48*nClusters;
    else %movie chunks that include primary alone patterns and primary+secondary patterns
        length = 50*nClusters;
    end
    Patterns = Pattern(i, 1:length);
    Times = TimeShift : Delay : TimeShift+Delay*(length-1);
    Chunk = NS_MovieChunkGenerationForExperiment(Times, nSamples, Patterns);
    MovieChunks = [MovieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

PatternOrder=Pattern;
MovieChunksFile=MovieChunks;

% from original version of function
% MovieChunks=[4];
% for i=1:4
%     l1=48*nClusters;
%     l2=58*nClusters;
%     if i<3
%         l=48*nClusters;
%     else
%         l=50*nClusters;
%     end
%     Patterns=Pattern(i,1:l);
%     Times=[TimeShift:Delay:TimeShift+Delay*(l-1)];
%     Chunk=NS_MovieChunkGenerationForExperiment(Times,nSamples,Patterns);
%     MovieChunks=[MovieChunks Chunk];
% end

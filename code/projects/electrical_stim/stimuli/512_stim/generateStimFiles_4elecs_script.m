% Script to generate 4 electrode stimulation patterns on the 512-electrode
% array

% L Grosberg 2014-07-15

% choose 1-3 electrodes (can be on edge of the array, but can't be 2nd nearest-neighbors)
% Requirements for clusters of electrodes: must be non-overlapping

%% electrode groups (pick 4 elecs), up to 3 clusters
% electrodeGroups = [164 155 156 148];%; 168 126 160 123]; % size(electrodeGroups,2) must be 4
% electrodeGroups = [164 155 156 148; 168 126 160 123]; % size(electrodeGroups,2) must be 4
% electrodeGroups  = [164 155 156 148; 78 70 75 67; 466 467 474 475]; % size(electrodeGroups,2) must be 4
% electrodeGroups  = [213 214 205 206; 264 272 260 268; 77 69 82 74]; % size(electrodeGroups,2) must be 4
% electrodeGroups  = [300 308 81 73; 45 37 42 34; 338 346 343 351]; % size(electrodeGroups,2) must be 4
% electrodeGroups  = [342 350 338 346; 465 466 473 474; 490 491 497 498]; % size(electrodeGroups,2) must be 4
electrodeGroups  = [3 6 11 473; 170 171 162 163; 461 462 468 461]; % size(electrodeGroups,2) must be 4

numClusters      = size(electrodeGroups,1); 
timeShiftInMs    = 0;        % offset the first stimulus from the beginning movie chunk
arraySpacing     = 60;       % 30 or 60 micron spacing
delayInMs        = 40;       %    
clusterDelayInMs = 10;       % clusterDelayInMs*numClusters < delayInMs
numberOfSamples  = 10000;    % 0.5s SET TO TRIGGER INTERVAL IN LABVIEW

if clusterDelayInMs*numClusters > delayInMs
    disp('Cluster stimulation will not complete within the desired rep rate for a single cluster'); 
    return; 
end
electrodeGroups  = electrodeGroups'; 
electrodes       = electrodeGroups(:); 
fprintf('Cluster 1:  %d %d %d %d\nCluster 2: %d %d %d %d\nCluster 3: %d %d %d %d\n',electrodes(1:4),electrodes(5:8),electrodes(9:12));
angleIncrement   = 2*pi/12;  % 30 degrees
% angleIncrement   = 2*pi/16; 
[Array, ~] = patternGen4D(angleIncrement);
nPatterns = size(Array,2); 

% Check to see if any patterns are repeated
if 0
for k = 1:size(Array,2)
    check = Array(:,k);
    m = zeros(1,size(Array,2));
    for i = 1:size(Array,2)
        if all(Array(:,i) == check)
            m(i) = 1;
        end
    end
    if size(find(m),2)>1
        disp(['k= ' num2str(k) ' i= ' num2str(size(find(m),2))]);
    end
end
end

% Generate movie files
%convert inputs to samples
Delay          = delayInMs*20;
TimeShift      = timeShiftInMs*20; % in sampling periods (50 microseconds)
clusterDelay   = clusterDelayInMs*20; 

if numClusters == 1
    Times          = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift; 
    nChunksTotal    = ceil(nPatterns/size(Times,2));
    movieChunks   = nChunksTotal;   % first value indicates number of chunks
    len          = size(Times,2);
    Pattern      = [1:1:nPatterns 1:(len*movieChunks-nPatterns)]; 
    Pattern      = reshape(Pattern, movieChunks, []);
elseif numClusters == 2
    Array         = blkdiag(Array,Array);
    nPatterns     = nPatterns*2; 
    
    Times1        = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift; 
    Times2        = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift + clusterDelay; 
    Times         = sort(cat(2,Times1,Times2)); 
    nChunksTotal  = ceil(nPatterns/size(Times,2));
    movieChunks   = nChunksTotal;   % first value indicates number of chunks
    len           = size(Times,2);
     
    Pattern       = 1:nPatterns; 
    a             = Pattern(:,1:nPatterns/2)';
    b             = Pattern(:,nPatterns/2+1:end)';
    Pattern       = reshape([a b]',size(a,1) + size(b,1),[])';
    Pattern       = [Pattern Pattern(1:(len*movieChunks-nPatterns))];
    Pattern       = reshape(Pattern, movieChunks, []);
elseif numClusters == 3
    Array         = blkdiag(Array,Array,Array);
    nPatterns     = nPatterns*3; 
    
    Times1        = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift; 
    Times2        = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift + clusterDelay;
    Times3        = (0:Delay:(numberOfSamples-Delay-TimeShift)) + TimeShift + clusterDelay*2; 
    
    Times         = sort(cat(2,Times1,Times2,Times3)); 
    
    nChunksTotal  = ceil(nPatterns/size(Times,2));
    movieChunks   = nChunksTotal;   % first value indicates number of chunks
    len           = size(Times,2);
    
    Pattern       = 1:nPatterns; 
    
    Pattern       = reshape(reshape(Pattern,[],numClusters)',1,[]); 
    Pattern       = [Pattern Pattern(1:(len*movieChunks-nPatterns))];
    Pattern       = reshape(Pattern, movieChunks, []);
else
    disp('undefined for more than 3 clusters of electrodes'); 
    return; 
end

% disp(['record for ' num2str(movieChunks*numberOfSamples/20000) ' s']); 
%%
for i =1:nChunksTotal
    Patterns = Pattern(i,1:len);
    Chunk=NS_MovieChunkGenerationForExperiment(Times, numberOfSamples, Patterns);
    movieChunks = [movieChunks Chunk]; %#ok<AGROW> %first value indicates number of chunks
end

movieChunksFile=movieChunks;

keyboard;
fid = fopen('4elec_3cluster_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('4elec_3cluster_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double'); 
fclose(fid);

fid = fopen('4elec_3cluster_movie','wb','ieee-le.l64');
fwrite(fid,movieChunksFile,'int32');
fclose(fid);
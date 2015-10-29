%% Safe zone combination tests
% Karthik's 'safe zone' calculations uses the single electrode stimulation
% to find the stimulation currents that activate an axon bundle. This
% script tests how neighboring pairs of electrodes around those threshold
% currents combine to give new 'safe zone' thresholds. 
% L Grosberg 9/24/2015
%% Define electrodes variable
saveFiles = 0; 
saveName = 'safeZone_2015-10-06-x';
% centerElectrodes = [272 107 79]; % Choose 3 for now
centerElectrodes = [376 45 304]; % Choose 3 for now

electrodes = [];
% Get the surrounding electrodes
for e = 1:length(centerElectrodes)
    currentCluster = getCluster512(centerElectrodes(e));
    electrodes = [electrodes currentCluster]; %#ok<AGROW>
end

% read the list of thresholds for electrode vector
% load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/safelevels.mat')
% safeZoneThresholds = safelevels(electrodes); 
% safeZoneThresholds = rand(size(electrodes))*4;  %Random vals for testing.
% thresholds = [1.71 1.61 1.31];
%    elec num         376    368   371   372   379   380   384    45    37    41    42    49    50    53   304   296   299   300   307   308   312
safeZoneThresholds = [1.71   3.7   1.91  1.1   1.31  1.31  .99    1.61  .99   .88   1.61  2.5   1.71  1.2  1.31  2.81  2.11  1.91  1.61  2.31  1.91]';

if length(unique(electrodes)) ~= length(electrodes)
    disp('Overlapping electrode choices!');
    return;
end
%% Define array variable
sz1 = safeZoneThresholds(1); 
sz2 = safeZoneThresholds(2); 
radius = norm([sz1 sz2],2); %Distance to corner of quadrant in space.
angle = atand(sz2/sz1); % Angle of the ray that passes through corner

sz1_start = sz1/2; % x-limit of the fake "quadrant"
sz2_start = sz2/2; % y-limit of the fake "quadrant"

% Get the angles of the extremal directions defined by fake "quadrant"
extremeangle1 = atand(sz2/sz1_start);
extremeangle2 = atand(sz2_start/sz1);

% Define the angles for the sweep (15 points between the extreme angles)
theta = linspace(extremeangle1, extremeangle2, 16);

% Define the x coordinates and y coordinates corresponding to these pts
e1 = (radius)*cosd(theta); %x coords
e2 = (radius)*sind(theta); %y coords
figure; scatter(e1,e2); 
line([sz1 sz1],get(gca,'YLim')); 
line(get(gca,'XLim'),[sz2 sz2]); 

%%
array = [];
for a = 1:size(theta,2)
    array2 = zeros(length(currentCluster),21);
    startPat = 1;
    for ii = 1:length(currentCluster)-1
        sz1 = safeZoneThresholds(ii); 
        sz2 = safeZoneThresholds(ii+1); 
        radius = norm([sz1 sz2],2);
        sizeOfEye = length(currentCluster)-ii;
        endPat = startPat + sizeOfEye - 1;
        array2(ii,startPat:endPat) = (radius)*cosd(theta(a));
        array2((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = (radius)*sind(theta(a))*eye(sizeOfEye);
        startPat = endPat + 1;
        array2
    end
    array = [array array2];
end
array = blkdiag(array,array,array);
figure; imagesc(array)
%% 2015-09-23 array variable
% 
% array = zeros(length(currentCluster),length(currentCluster)*(length(currentCluster)+1)/2 - length(currentCluster)); 
% startPat = 1;
% for ii = 1:length(currentCluster)-1
%     sizeOfEye = length(currentCluster)-ii;
%     endPat = startPat + sizeOfEye - 1; 
%     array(ii,startPat:endPat) = 1;
%     array((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = eye(sizeOfEye);
%     startPat = endPat + 1; 
% end
% array2 = zeros(length(currentCluster),21); 
% startPat = 1;
% for ii = 1:length(currentCluster)-1
%     sizeOfEye = length(currentCluster)-ii;
%     endPat = startPat + sizeOfEye - 1; 
%     array2(ii,startPat:endPat) = 1;
%     array2((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = 1.1*eye(sizeOfEye);
%     startPat = endPat + 1; 
% end
% array = [array array2];
% 
% array2 = zeros(length(currentCluster),21); 
% startPat = 1;
% for ii = 1:length(currentCluster)-1
%     sizeOfEye = length(currentCluster)-ii;
%     endPat = startPat + sizeOfEye - 1; 
%     array2(ii,startPat:endPat) = 1;
%     array2((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = 1.1*1.1*eye(sizeOfEye);
%     startPat = endPat + 1; 
% end
% 
% array = [array array2];
% 
% array2 = zeros(length(currentCluster),21); 
% startPat = 1;
% for ii = 1:length(currentCluster)-1
%     sizeOfEye = length(currentCluster)-ii;
%     endPat = startPat + sizeOfEye - 1; 
%     array2(ii,startPat:endPat) = 1;
%     array2((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = eye(sizeOfEye)/1.1;
%     startPat = endPat + 1; 
% end
% 
% array = [array array2];
% 
% array2 = zeros(length(currentCluster),21); 
% startPat = 1;
% endPat = (length(currentCluster)-1); 
% for ii = 1:length(currentCluster)-1
%     sizeOfEye = length(currentCluster)-ii;
%     endPat = startPat + sizeOfEye - 1; 
%     array2(ii,startPat:endPat) = 1;
%     array2((ii+1):length(currentCluster), startPat:1:(startPat+sizeOfEye)-1) = eye(sizeOfEye)/1.1/1.1;
%     startPat = endPat + 1; 
% end
% 
% 
% 
% array = [array array2];
% array = blkdiag(array,array,array);
% figure; subplot(2,1,1); imagesc(array)
% arrayScaled = array.*repmat(safeZoneThresholds,1,size(array,2));
% subplot(2,1,2);  imagesc(arrayScaled);
% ylabel('electrode'); xlabel('pattern number'); 
%%
TimeShiftInMs=0;
InterPulseLatencyInMs=25;
TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000; %trigger interval = 0.5 s
Patterns = reshape(1:size(array,2),[],length(centerElectrodes))'; 
Patterns = Patterns(:)'; 
Times=TimeShift:InterPulseLatency:TimeShift+(length(Patterns)-1)*InterPulseLatency;

nMovieChunks = ceil(Times(end)/NumberOfSamples);
allPatterns=[];
patternchunklength = ceil(length(Patterns)/nMovieChunks); 
Chunks=[];
for i=1:nMovieChunks; 
    Start=(i-1)*patternchunklength;
    try
        PatternsForMovie=Patterns(Start+1:Start+patternchunklength);
    catch
        PatternsForMovie=Patterns(Start+1:end);
    end
    mChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    allPatterns=[allPatterns mChunk];
end
MovieChunksFile=[nMovieChunks allPatterns];
%%
% Saving
if saveFiles
    fid = fopen([saveName '_mv'],'wb');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_el'],'wb');
    fwrite(fid,electrodes,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_pt'],'wb','ieee-le.l64');
    fwrite(fid,arrayScaled,'double');
    fclose(fid);
end


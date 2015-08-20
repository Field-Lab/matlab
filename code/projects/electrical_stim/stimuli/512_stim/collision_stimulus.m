% Collision stimulus used with Karthik during 2015-04 experiments
% L Grosberg April 2015

%% Manually enter electrodes that will be used for simultaneous stimulation
% Select sets on either side of an axon bundle (defined by Karthik's graph
% code) that will activate that bundle. Each electrode from set 1 will be
% paired with each from set 2 to simulataneously stimulate. All electrodes
% in each group (sets 1 and 2) will be activated simulateously in addition.
% 
% Each group defined stimulates a different bundle.

group(1).set1 = [412 413 414 415]; 
group(1).set2 = [96 88 80 72];  
group(2).set1 = 401:404;  
group(2).set2 = 132:135; 
group(3).set1 = [261 269 277 285]; 
group(3).set2 = 193:196; 
saveFiles = 0; % Set to 1 to save the stimulation files, 0 for testing
saveName = 'collision'; %Descriptive save name for the el, mv, and pt files

% Concatenate electrode sets to 
elecs = []; 
arraysize = size(group,2); % Initialize at the number of groups to include all electrodes in each group stimulating at once
for j = 1:size(group,2)
    elecs = [elecs group(j).set1 group(j).set2];  %#ok<AGROW>
    arraysize = arraysize+size(group(j).set1,2)*size(group(j).set2,2); 
end
array = zeros(size(elecs,2),arraysize); 
colBegin = 1; 
rowBegin = 0;
for j = 1:size(group,2)
    array(rowBegin+(1:(size(group(j).set1,2)+size(group(j).set2,2))),colBegin) = 1; 
%     rowBegin+(1:(size(group(j).set1,2)+size(group(j).set2,2)))
    for k = 1:size(group(j).set1,2)
        array(rowBegin + k, colBegin + (1:size(group(j).set2,2))) = 1; 
        array(rowBegin + size(group(j).set1,2)+(1:size(group(j).set2,2)) , colBegin + (1:size(group(j).set2,2))) = eye(size(group(j).set2,2)); 
        colBegin = colBegin+size(group(j).set2,2); 
    end
    colBegin = colBegin + 1; 
    rowBegin = rowBegin + size(group(j).set1,2)+size(group(j).set2,2); 
end
figure; imagesc(array); 
xlabel('pattern number'); ylabel('elecs used: see elecs variable for values'); 

TimeShiftInMs=5;
InterPulseLatencyInMs=25;
TimeShift=TimeShiftInMs*20;
InterPulseLatency=InterPulseLatencyInMs*20;

NumberOfSamples=10000; %trigger interval = 0.5 s
Patterns = reshape(1:size(array,2),[],size(group,2))'; 
Patterns = Patterns(:)'; 
Times=[TimeShift:InterPulseLatency:TimeShift+(length(Patterns)-1)*InterPulseLatency];

AllPatterns=[];
patternchunklength = size(array,2)/size(group,2); 
Chunks=[];
for i=1:ceil(Times(end)/NumberOfSamples) 
    Start=(i-1)*patternchunklength;
    PatternsForMovie=Patterns(Start+1:Start+patternchunklength);
    mChunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,PatternsForMovie);
    AllPatterns=[AllPatterns mChunk];
end
MovieChunksFile=[ceil(Times(end)/NumberOfSamples) AllPatterns];

% Saving
if saveFiles
    fid = fopen([saveName '_mv'],'wb');
    fwrite(fid,MovieChunksFile,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_el'],'wb');
    fwrite(fid,elecs,'int32');
    fclose(fid);
    
    fid = fopen([saveName '_pt'],'wb','ieee-le.l64');
    fwrite(fid,array,'double');
    fclose(fid);
end

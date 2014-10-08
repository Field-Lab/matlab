% Code to produce stimuli on four electrodes with random current
% assignments on each electrode in the pattern. 


TimeShiftInMs = 0;
DelayInMs     = 7.5; %time between subsequent pulses on different electrodes

NumberOfSamples = 40000;     % 2s SET TO TRIGGER INTERVAL IN LABVIEW
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


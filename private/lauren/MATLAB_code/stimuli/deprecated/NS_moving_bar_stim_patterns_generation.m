
NumberOfSamples=40000;
%electrodes=[1:8 10:24 26:56 58:64];
electrodes=[1:64];
Array=eye(64);
pathToNeuronsFile = '/Volumes/War/Analysis/Lauren/2008-04-22-2/data001-map/data001-map.neurons';

% ID numbers of chosen cells as given in VISION analysis (specified in neurons file)
cellIDs = [95 97 260 487];
nPatterns = length(cellIDs); %to generate dummy pattern numbers

%beginning and end of time window withing neurons file to extract activity pattern from
timeWindowStart = 1;
timeWindowEnd = 3;


% electrodeOrder1 = LaurenelectrodeOrderGenerator(0);
% electrodeOrder2 = LaurenelectrodeOrderGenerator(0);

% dummy pattern numbers, for use until we have mapping between stim pattern and cell ID
patternNumbers=1:nPatterns;

% retrieves spike times for specified neurons within specified time window (in samples, relative to
% timeWindowStart)
spikeTimesInWindow = getNeuronSpikeTimes(pathToNeuronsFile, cellIDs, timeWindowStart, timeWindowEnd);

[Times Patterns] = combineNeuronSpikeTimes(spikeTimesInWindow, patternNumbers);

Chunk=NS_MovieChunkGenerationForExperiment(Times, NumberOfSamples, Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd /Applications/MATLAB74/work/Lauren/stimuli/stimulus_files/; 
fid = fopen('movingbar_electrodes','wb')
fwrite(fid,electrodes,'integer*4');
fclose(fid);

fid = fopen('movingbar_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('movingbar_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);
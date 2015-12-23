%2011-08-31 removed 2ms rounding version of bar (1 ms is sufficient)

clear all

% pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2012-09-27-0/data003/data003.neurons';
% cellIDs =        [2794 2342 3293 1966];
% cellElectrodes = [14  31  39   47];
% stimAmps =       [1.7 1.6 1.25 1.3];


pathToNeuronsFile = '/Volumes/Analysis/2011-06-24-5-test/data001/data001.neurons';
cellIDs =        [1 32 183 333 381 783];
cellElectrodes = [1  2  14  23  26  50];
stimAmps =       [0.88 1.6 0.68 0.88 0.98 1.1];

%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data002-from-data000/data002-from-data000.neurons';

%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2012-01-27-2/data005/data005.neurons';
%cellIDs =          [198  334  904  483  724  621  755  812   695];
%cellElectrodes =   [18   23   1    28   38   45   51   52    41];
%stimAmps =         [5.3  3.1  2.8  2    3    4.8  5.4  4.4   2.3];


%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data001-from-data000/data001-from-data000.neurons';
% cellIDs =        [168 452 572  691  769  949];
% cellElectrodes = [14  31  39   47   50   63];
% %stimAmps = [0.65 1.7 2.1 1.6 3 0.8];
% stimAmps =       [1.7 1.6 1.25 1.3  1.7  1];


% pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data001-from-data000/data001-from-data000.neurons';
% cellIDs =          [2   33   183   393   604   738];
% cellElectrodes =   [60  5    14    21    26    48];
% stimAmps =         [0.9 1.9  0.65   0.7  0.8   0.73];

%pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-07-05-0/data003/data003.neurons';
%cellIDs =          [77   257  453  557  542  797];
% cellElectrodes =   [1    23   31   33   40   49];
% stimAmps =         [1.4  1.5  0.6  1.1  1.7  2.2];


% pathToNeuronsFile = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data001/data001.neurons';
% cellIDs =          [213  200  348  467  543  663  783  901];
% % cellElectrodes =   [4    16   24   28   34   45   54   61];
% % stimAmps =         [2.2  2.0  1.3  0.8  4.1  1.3  2.8  1.4];


nSamples = 20000; %make sure this is enough time for moving bar responses


stimulusType = 'MB';  %'NS'; % 'MB' (moving bar) or 'NS' (natural scenes)

%%%%%%%%%%%% parameters for visual stimulus %%%%%%%%%%%%%%%%%%%%%
if strcmpi(stimulusType, 'MB') %LG nts checking for trigger structures
    barColor = 'white'; % 'white' or 'black'
elseif strcmpi(stimulusType, 'NS')
    nScenes = 5;
    sceneLength = 1; %in seconds
    sceneIndex = 1;
else
    error('stimulus type not recognized')
end

%should match what is set on lisp code
nReps = 100; % number of repetitations of chosen stimulus (e.g. one color of moving bar or one natural scene)


%%%%%%%%%%%% parameters for natural scenes stimulus %%%%%%%%%%%%%%%%%%%%%


% Chooses a 'representative' response to the visual stimulus
nToChoose = 2; %number of repetitions to pick out for replicating
choiceStart = 1; %in the ordering of least-to-most mean distance to other spike trains, start with this one (in most cases, = 1)

% cellIDs =        [36   48   247  393  531  694];
% cellElectrodes = [2    4    17   21   36   47];
% stimAmps =       [1.1  0.68 .35  1.1  .55  0.9];

% adds individual electrode stims to overal stimulus file to compare
% threshlds for isolated pulses vs pulses in a spatiotemporal seq
includeIndividualStim = true; %whether or not to add a final movie chunk that applies each single-electrode stimuli individually

%startPadding = 0.45; %how long after start trigger response is visible (in seconds)
%endPadding = 0.4155; %how long before end trigger reponses seems to end (in seconds)

% Trial and error basis to figure out when bar actually passing over array
startPadding = 0.15; % in seconds
endPadding = 0.05;

%used in spike-distance metric talk to EJ
costParam = 10; %in ms, the distance that a spike moves that is equivalent (in cost) to deleting and inserting a spike

normType = 'both'; %type of normalization (can be 'none', 'origin', 'destination', 'both')

%%

if strcmpi(stimulusType, 'MB')
    [spikeTimesAll intervals] = loadMovingBarResponses(pathToNeuronsFile, cellIDs, barColor, 'trialPadding', [startPadding endPadding]);
else %natural scenes
    [spikeTimesAll intervals] = loadNaturalSceneResponses(pathToNeuronsFile, cellIDs, nScenes, sceneLength, sceneIndex, 'trialPadding', [startPadding endPadding]);
end

if length(spikeTimesAll) ~= nReps
    warndlg('unexpected number of stimulus repetitions--check that nReps value is correct')
end

windowLength = intervals{1}(2) - intervals{1}(1);


[spikeTimesToUse order] = spike_train_median_calculator(spikeTimesAll, nToChoose, choiceStart, costParam, normType, windowLength, 'cellIDs', cellIDs);


%% plots chosen spike trains 

if 1
    nCells = length(cellIDs);
    figure
    for i = 1:nCells
        axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
        hold on
        for k = 1:nToChoose
            for j = 1:length(spikeTimesToUse{k}{i})
                plot([spikeTimesToUse{k}{i}(j) spikeTimesToUse{k}{i}(j)], [k-1 k], 'k-', 'LineWidth', 1)
            end
        end
        hold off

        set(gca, 'yLim', [-5 nToChoose+5], 'xlim', [0 nSamples/20000])
        if i == nCells
            xlabel('time (s)')
            set(gca, 'ylim', [0 nToChoose])
        else
            set(gca, 'xtick', [], 'ylim', [0 nToChoose])
        end
        ylabel(['elec ' num2str(cellElectrodes(i)), 10, 'amp ' num2str(stimAmps(i))])
    end
end

%% generate stimulus files

electrodes = cellElectrodes;
Array = diag(stimAmps);


%convert to samples
for ii = 1:nToChoose
    for jj = 1:length(spikeTimesToUse{ii})
        spikeTimesToUse{ii}{jj} = spikeTimesToUse{ii}{jj}*20000;
    end
end


MovieChunksFile = [];

% concatenate spike times for all cells and round to nearest 1 or 2 ms
for ii = 1:nToChoose %loops through chosen reps
    %spikeTimesOrigCat = [];
    spikeTimes1Cat = [];
    %spikeTimes2Cat = [];
    %patternsOrig = [];
    patterns1 = [];
    %patterns2 = [];
    for jj = 1:length(spikeTimesToUse{ii}) %loops through cells
        %spikeTimesOrig = floor(spikeTimesToUse{ii}{jj}); %rounded to nearest sample
        spikeTimes1 = 20*floor(spikeTimesToUse{ii}{jj}/20); %rounded to nearest 1 ms
        %spikeTimes2 = 40*floor(spikeTimesToUse{ii}{jj}/40); %rounded to nearest 2 ms

        %make sure 2 spikes haven't been rounded to the same value, and remove one if they did
        if any(diff(spikeTimes1)==0)
            spikeTimes1(find(diff(spikeTimes1)==0)) = [];
            warning('2 or more spikes rounded to the same 1 ms')
        end
%         if any(diff(spikeTimes2)==0)
%             spikeTimes2(find(diff(spikeTimes2)==0)) = [];
%             warning('2 or more spikes rounded to the same 2 ms')
%         end

        %spikeTimesOrigCat = [spikeTimesOrigCat; spikeTimesOrig];
        spikeTimes1Cat =    [spikeTimes1Cat;    spikeTimes1]; %#ok<AGROW>
        %spikeTimes2Cat =    [spikeTimes2Cat;    spikeTimes2]; %#ok<AGROW>
        %patternsOrig = [patternsOrig; jj*ones(length(spikeTimesOrig),1)]; %#ok<AGROW>
        patterns1 =    [patterns1;    jj*ones(length(spikeTimes1),1)]; %#ok<AGROW>
        %patterns2 =    [patterns2;    jj*ones(length(spikeTimes2),1)];
    end
    
    %make sure all spike times are within nSamples
    %if any(spikeTimesOrigCat >= nSamples) || any(spikeTimes1Cat >= nSamples) || any(spikeTimes2Cat >= nSamples)
    %if any(spikeTimesOrigCat >= nSamples) || any(spikeTimes1Cat >= nSamples)
    if any(spikeTimes1Cat >= nSamples)
        errordlg('nSamples is not large enough for entire moving bar response')
    end

    % sorts spikeTimesCat and patterns by time
    %[spikeTimesOrigCat iSortOrig] = sort(spikeTimesOrigCat);
    [spikeTimes1Cat    iSort1]    = sort(spikeTimes1Cat);
    %[spikeTimes2Cat    iSort2]    = sort(spikeTimes2Cat);
    %patternsOrig = patternsOrig(iSortOrig);
    patterns1    = patterns1(iSort1);
    %patterns2    = patterns2(iSort2);
    
    %ChunkOrig = NS_MovieChunkGenerationForExperiment(spikeTimesOrigCat, nSamples, patternsOrig);
    Chunk1 = NS_MovieChunkGenerationForExperiment(spikeTimes1Cat, nSamples, patterns1);
    %Chunk2 = NS_MovieChunkGenerationForExperiment(spikeTimes2Cat, nSamples, patterns2);
    
    %MovieChunksFile=[MovieChunksFile ChunkOrig Chunk1 Chunk2];
    %MovieChunksFile=[MovieChunksFile ChunkOrig Chunk1];
    MovieChunksFile = [MovieChunksFile Chunk1];
end

if includeIndividualStim
    patterns = 1:length(spikeTimesToUse{ii});
    times = 0:50*20:50*20*(length(spikeTimesToUse{ii})); %apply stimuli every 50 ms
    chunk = NS_MovieChunkGenerationForExperiment(times, nSamples, patterns);
    MovieChunksFile = [MovieChunksFile chunk];
    %MovieChunksFile = [3*nToChoose+1 MovieChunksFile];
    %MovieChunksFile = [2*nToChoose+1 MovieChunksFile];
    MovieChunksFile = [nToChoose+1 MovieChunksFile];
else
    %MovieChunksFile = [3*nToChoose MovieChunksFile];
    %MovieChunksFile = [2*nToChoose MovieChunksFile];
    MovieChunksFile = [nToChoose+1 MovieChunksFile];
end

keyboard

%% save out files

fid = fopen('moving_bar_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('moving_bar_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('moving_bar_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 
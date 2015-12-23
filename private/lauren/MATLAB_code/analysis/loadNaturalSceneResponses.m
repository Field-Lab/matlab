function [spikeTimesAll intervals] = loadNaturalSceneResponses(pathToNeuronsFiles, cellIDs, nScenes, sceneLength, sceneIndex, varargin)
%   extracts spike times for set of neurons in response to 'standard' moving bar
%   stimulus: alternating black and white bar (white first) that sends triggers at
%   beginning and end of each bar repetition
%
%   standard moving bar takes between 1s and 1.5s per repetition, so
%   warning is displayed if time between triggers is < 1s
%
%
% arguments:
%   pathToNeuronsFile: self explanatory
%   nScenes: number of different scenes that are applied in the sequence
%   sceneLength: length of one scene (in seconds)
%   sceneIndex: index of the scene to extract spikes for
%   cellIDs: vector of neuron IDs
%   startPadding: time between start of natural scene and beginning of extracted time window
%   endPadding: time between end of extracted time window and end trigger
%
% returns:
%   spikeTimesAll: cell array of spike times -- spikeTimesAll{movingBarRepIndex}{cellIndex}(spikeIndex)
%   intervals: cell array of time window boundaries -- intervals{movingBarRepIndex}(windowStart windowEnd)
%   

%% parse inputs

p = inputParser;

p.addRequired('pathToNeuronsFile', @ischar)
p.addRequired('nScenes', @isnumeric)
p.addRequired('sceneLength', @isnumeric)
p.addRequired('sceneIndex', @isnumeric)
p.addRequired('cellIDs', @isnumeric)

p.addParamValue('trialPadding', [0 0], @isnumeric)

p.parse(pathToNeuronsFiles, nScenes, sceneLength, sceneIndex, cellIDs, varargin{:})

startPadding = p.Results.trialPadding(1);
endPadding = p.Results.trialPadding(2);

%% get triggers to determine start/end times of moving bar repetitions


datarun.names.rrs_neurons_path = pathToNeuronsFile;
datarun = load_neurons(datarun);

triggers = datarun.triggers;

%determine expected number of triggers per repetition
repLength = sceneLength*nScenes;
trigPerRep = floor(repLength/0.8333)+1;

trigInts = diff(triggers);
trigIntsDuringStim = trigInts;
trigIntsDuringStim(trigPerRep:trigPerRep:end) = []; %remove intervals that correspond to time between repetitions
if any(trigIntsDuringStim(1) < 0.83) || any(trigIntsDuringStim(1) > 0.84)
    error('unexpected trigger timing: triggers should be every 100 frames (~833 ms) while natural scenes are applied')
end


tScene = sceneLength*(sceneIndex-1); %when scene starts relative to trigger marking start of repetition

repStarts = triggers(1:trigPerRep:end);
sceneStarts = repStarts + tScene;
sceneEnds = sceneStarts + sceneLength;

nReps = length(sceneStarts);


nCells = length(cellIDs);


%load all spike times from natural scene repetitions
spikeTimesAll = cell(nReps, 1);
intervals = cell(nReps, 1);

for i = 1:nReps
    spikeTimesAll{i} = cell(nCells, 1);
    intervals{i} = [sceneStarts(i)+startPadding sceneEnds(i)-endPadding];
    for j = 1:nCells
        spikeTimesTemp = datarun.spikes{get_cell_indices(datarun, cellIDs(j))};
        % retrieves spike times within window of time
        spikeTimesAll{i}{j} = spikeTimesTemp(spikeTimesTemp > intervals{i}(1) & spikeTimesTemp < intervals{i}(2));
        spikeTimesAll{i}{j} = spikeTimesAll{i}{j} - intervals{i}(1); % sets start of time window to time = 0
    end
end

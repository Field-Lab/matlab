function [spikeTimesAll intervals] = loadMovingBarResponses(pathToNeuronsFile, cellIDs, barColor, varargin)
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
%   /snle/lab/Experiments/Array/Analysis/2011-06-24-5/data006/data006.neurons
%   cellIDs: vector of neuron IDs
%   barColor: either 'white' or 'black'
%   startPadding: time between start trigger and beginning of extracted time window
%   endPadding: time between end of extracted time window and end trigger
%
% returns:
%   spikeTimesAll: cell array of spike times -- spikeTimesAll{movingBarRepIndex}{cellIndex}(spikeIndex)
%   intervals: cell array of time window boundaries -- intervals{movingBarRepIndex}(windowStart windowEnd)
%   

%% parse inputs

p = inputParser;

p.addRequired('pathToNeuronFile', @ischar)
p.addRequired('cellIDs', @isnumeric)
p.addRequired('barColor', @ischar)

p.addParamValue('trialPadding', [0 0], @isnumeric)
p.addParamValue('expectedLengthLims', [1 2], @isnumeric)

p.parse(pathToNeuronsFile, cellIDs, barColor, varargin{:})


startPadding = p.Results.trialPadding(1);
endPadding = p.Results.trialPadding(2);
intLims = p.Results.expectedLengthLims;

%% get triggers to determine start/end times of moving bar repetitions


datarun.names.rrs_neurons_path = pathToNeuronsFile;
datarun = load_neurons(datarun);

triggers = datarun.triggers;

if strcmpi(barColor, 'white')
    triggersEnd = triggers(2:4:end); %end of moving bar
    nBarReps = length(triggersEnd);
    triggersStart = triggers(1:4:nBarReps*4-3); %start of moving bar
elseif strcmpi(barColor, 'black')
    triggersEnd = triggers(4:4:end); %end of moving bar
    nBarReps = length(triggersEnd);
    triggersStart = triggers(3:4:nBarReps*4-1); %start of moving bar
elseif strcmpi(barColor, 'all')
    triggersEnd = triggers(2:2:end); %end of moving bar
    nBarReps = length(triggersEnd);
    triggersStart = triggers(1:2:nBarReps*2-1); %start of moving bar
else
    errordlg('invalid bar color specified')
end

if min(triggersEnd - triggersStart) < intLims(1) || max(triggersEnd - triggersStart) > intLims(2)
    warndlg(['moving bar takes less than ' num2str(intLims(1)) ' second or more than ' num2str(intLims(2)) '--check that triggers are identified correctly'])
end

nCells = length(cellIDs);


%load all spike times from moving bar repetitions
spikeTimesAll = cell(nBarReps, 1);
intervals = cell(nBarReps, 1);

for i = 1:nBarReps
    spikeTimesAll{i} = cell(nCells, 1);
    intervals{i} = [triggersStart(i)+startPadding triggersEnd(i)-endPadding];
    for j = 1:nCells
        spikeTimesTemp = datarun.spikes{get_cell_indices(datarun, cellIDs(j))};
        % retrieves spike times within window of time
        spikeTimesAll{i}{j} = spikeTimesTemp(spikeTimesTemp > intervals{i}(1) & spikeTimesTemp < intervals{i}(2));
        spikeTimesAll{i}{j} = spikeTimesAll{i}{j} - intervals{i}(1); % sets start of time window to time = 0
    end
end

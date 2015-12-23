function break_up_neurons(neuronsPath, timeRange, varargin)
% BREAK_UP_NEURONS  Break up a standard neurons file into multiple pieces
% specified by time intervals in the timeRange vector. These neurons files
% will be saved to neuronsPath by default, although an alternative path can
% be specified.
%
% usage: break_up_neurons(neuronsPath, timeRange, <params>)
%
% arguments: neuronsPath - path to neurons file to be broken up
%            ex: '/Analysis/Machado/2005-04-26-0/data000/data000.neurons'
%            
%            timeRange - cell array of time intervals (in seconds)
%            ex:   range{i} = [begin end];
%                  range{1} = [60  16*60]; range{2} = [16*60 31*60];
%                        NOTE: time intervals that start or extend beyond
%                        the length of the neurons file will be ignored.
%
%                              -OR-
%
%            timeRange - length of each interval (in seconds)
%            ex: timeRange = 60*10 will break the whole run into 10 minute
%            chunks and will figure out the bounds for each interval
%            automatically.
%            
%            params - struct or list of optional parameters (see below)
%
% output:    writes neurons files specified by timeRange to outputPath
%
%
% optional params, their default values, and what they specify:
%
% verbose           false         show output text
% outputPath        neuronsPath   save neurons file to this location
% samplingRate      20000         if this changes, vision will break anyway
% createSubfolders  true          create subfolders for each neurons file
%
%
% 4/23/08 tamachado
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('outputPath', neuronsPath);
p.addParamValue('samplingRate', 20000);
p.addParamValue('createSubfolders', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% set up timer
if params.verbose
    fprintf('\n Breaking up neurons file...');
    start_time = clock; % note when it started
end

%%%%%%%%%%%%%%%%%%%%%%
% Load up neurons file
%%%%%%%%%%%%%%%%%%%%%%

% open up the neurons file
nf = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronsPath);

% get the list of cell ids
cellIds = nf.getFullIDList;

% check if input range is a cell or if we need to make it a cell
if ~iscell(timeRange)
    nSeconds = nf.getNumberOfSamples/params.samplingRate;
    range = cell(ceil(nSeconds/timeRange),1);
    interval = floor(nf.getNumberOfSamples/length(range));
    
    for rng  = 0:length(range)-1
        range{rng+1} = [interval*rng interval*(rng+1)]+1;
    end
    timeRange = range;
end

% for each cell, parcel it up into n cells where
% n is the number of output neurons files we want
spikeTimes = cell(length(cellIds),length(timeRange));
electrodes = zeros(length(cellIds),1);
emptySection = ones(length(timeRange),1);
    
    
for id = 1:length(cellIds)
    % store all spike times for current cell
    st = nf.getSpikeTimes(cellIds(id));
    
    % store electrode where id was found
    electrodes(id) = nf.getElectrode(cellIds(id));
    
    % save out the specified spike times
    for rng = 1:length(timeRange)
        r = timeRange{rng};
        if r(1) > r(2), continue; end
         spikeTimes{id,rng} = st(intersect(find(st >= r(1)),find(st <= r(2))));
        
        % keep track of which segments are empty
        if sum(spikeTimes{id,rng}) > 0
            emptySection(rng) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out neurons files
%%%%%%%%%%%%%%%%%%%%%%%%%

% remove the file name from the output path (if there is one)
dotIndex = strfind(neuronsPath,'/');
if ~isempty(dotIndex)
    neuronsPath(dotIndex(end):end) = [];
end

% put all pieces in a subfolder called pieces
if params.createSubfolders
    neuronsPath = [neuronsPath '/pieces'];
    unix(['mkdir ' neuronsPath]);
end

% create each neuron file
for rng = 1:length(timeRange)
    % name this output file
    fn = sprintf('part%d',rng);
    
    % if this segment is empty, skip it
    if emptySection(rng)
        continue;
    end
    
    % create a subfolder for each neuron file (if desired)
    if params.createSubfolders
        newPath = [neuronsPath '/' fn];
        unix(['mkdir ' newPath]);
    else
        newPath = neuronsPath;
    end
    
    % write the empty file to disk
    newNf = edu.ucsc.neurobiology.vision.matlab.ReadNeurons([newPath '/' fn '.neurons'], nf);
    
    % populate it with each cell
    for id = 1:length(cellIds)
        % if a given cell is empty on this segment, don't add it
        % vision crashes if there are 0 or 1 spikes in a given cluster
        if length(spikeTimes{id,rng}) < 2
            continue;
        end
        
        % add the current cell to this segment
        newNf.addNeuron(electrodes(id), spikeTimes{id,rng}, cellIds(id), length(spikeTimes{id,rng}));
    end
end

% create a globals file for each neurons file
if params.createSubfolders
    makeGlobals = edu.ucsc.neurobiology.vision.convert.AddGlobalsFiles;
    makeGlobals.topLevel(neuronsPath,true);
end

% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end


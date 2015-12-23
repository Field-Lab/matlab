function [spikes, extras] = load_rrs_neurons(full_path, neurons)
% LOAD_RRS_NEURONS    loads up RRS .neurons file
%
% usage:   [spikes] = load_rrs_neurons(neurons_file);
%          [spikes, extras] = load_rrs_neurons(...);
%          [...] = load_rrs_neurons(neurons_file, neurons);
%
% arguments:  neurons_file - full path of .neurons file
%                  neurons - [optional]
%
%                            If vector of neurons IDs, then only
%                            extract particular spike times.
%
%                            If argument not specified, or is 'all', 
%                            then extract all neurons. 
%                            
%                            If empty vector, no spike times extracted
%                            but all cell ids and channels returned.
%
% outputs:    spikes       - Nx1 cell array where N cells loaded.
%                            Each element is vector of spike times (s).
%
%             extras       - structure with lots of optional stuff
%                             | duration  -- total duration (s)
%                             | triggers  -- Nx1 vector of trigger times
%                             | channels  -- Nx1 vector of electrodes
%                             | cell_ids  -- Nx1 vector of cell ID's
%
% Reads .neurons file used by SNL-E and SCIPP software. Routines
% were largely based off of EJC's LISP code "read-rrs-cell-file"
% for reading in cell files.
%
% examples:  
%     spikes = load_rrs_neurons('/Data/test/data000.neurons','all');
%     [spikes, extras] = load_rrs_neurons('/Data/test/data000.neurons');
%     [spikes, extras] = load_rrs_neurons('/Data/data000.neurons',[67]);
%
% see also:  LOAD_NEURONS, LOAD_MMM_NEURONS
%
%  2006-08-19 shlens@salk.edu  Minor bug fixes.
%  2005-06-21 shlens@salk.edu
%

%% 1. Definitions %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data types
long_type = 'int32';
unsigned_long_type = 'uint32';
float_type = 'float64';
double_type = 'double';
long_size = 4;

% constants
unused_slot_tag = -(2^31);
end_of_header = 152;
blank_value = -2;

%% 2. Open File   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~((nargin == 1) || (nargin == 2))
  error('load_rrs_neurons: incorrect number of arguments');
end
  
% open file in appropriate format
fid = fopen(full_path,'r','ieee-be');

if (fid == -1)
  error('load_rrs_neurons: did not find file %s',full_path);
end


%% 3. Read Header %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in header information
header = fread(fid,4,long_type);

file_version = header(1);        % file version
header_slots = header(2);        % headerCapacity (max # of cells)
nTicks= header(3);               % nSamples (in ticks)
sampling_frequency = header(4);  % Hz

% check sampling_frequency
if sampling_frequency ~= 20000
  error('load_rrs_neurons: incorrect sample frequency')
end

% file creation info (not used)
time_hi = fread(fid,1,unsigned_long_type);
time_lo = fread(fid,1,unsigned_long_type);


%% 3. Determine Version Params %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default in case parameters were not set
contamination = 0;
min_spikes = 0;
remove_duplicates = 0;

switch file_version
 case 32
  %% --- Very old Obvius version
  %% --  This is the version used by all Java code
  number_of_records = 1; 
  spike_time_type = long_type;
 
  % extra parameters (2005-11-05)
  contamination     = fread(fid,1,double_type);
  min_spikes        = fread(fid,1,long_type);
  remove_duplicates = fread(fid,1,long_type) == 1;
  
  % check if blank values
  if (min_spikes == blank_value) min_spikes = 0; end;
  if (contamination == blank_value) contamination = 0; end;
  
 case 33
  %% --- Current Obvius version (unknown?)
  number_of_records = fread(fid,1,long_type);
  spike_time_type = long_type;
 case 100
  %% -- Dumitru's version from Vision used in Manual sorting
  number_of_records = 1;
  spike_time_type = float_type;
 otherwise
  error('unknown version')
end


% skip through end of blank header
fseek(fid, end_of_header, -1);


%% 4. Cell ID Header Information %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make something [4 x header_slots]
buffer = zeros(header_slots,4);
buffer = fread(fid, [4 header_slots], long_type);

% extract useful info
cell_ids = buffer(1,:)';
channels = buffer(2,:)';

% total number of cells (subtract 1 for lisp version)
if isempty(find(channels == unused_slot_tag))
  % version 33 uses this method for finding number of cells
  number_of_cells = length(cell_ids);
else
  % version 32, 100 uses this method for finding number of cells
  number_of_cells = min(find(channels == unused_slot_tag)) - 1;
end

% remove unused slots (vast majority)
cell_ids = cell_ids(1:number_of_cells);
channels = channels(1:number_of_cells);

% some channels are read as negative numbers (why??)
channels = abs(channels);

% subtract out duplicate cells
% num_cells = number_of_cells - length(find(channels < 0));

% check first cell is trigger
if channels(1) ~= 0
  error('load_rrs_neurons: trigger not found')
end

% delete buffer
clear buffer

% update user
fprintf(1,['Examining ' num2str(number_of_cells) ' cells (RRS v.' ...
	   num2str(file_version) ') ... ']);


%% 5. Determine neurons to extract %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% neuron ID for triggers
trigger_id = cell_ids(1);

% determine user choice
if (nargin == 1) || ischar(neurons)
  % default: load all neurons
  neurons = cell_ids;
else
  % always load the triggers
  neurons = reshape(neurons, length(neurons), 1);
  neurons = [trigger_id; neurons];
end  

% find if any neurons requested do not exist
intersection = intersect(cell_ids, neurons);
bad_neurons = setdiff(neurons, intersection);
if ~isempty(bad_neurons)
  error(['load_rrs_neurons: failed to find neurons ['  ...
	 num2str(bad_neurons') ']']);
end

% ensure that neurons are unique
if length(neurons) ~= length(unique(neurons))
  error(['load_rrs_neurons: duplicate copy of neuron requested']);
end

%% 6. Allocate Memory  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes = cell(length(neurons),1);
electrodes = zeros(length(neurons),1);
spike_counts = zeros(number_of_cells,1);


%% 7. Load Spike Times  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through the records, then each cell
for i=1:number_of_records

  % keep track of number of neurons extracted
  neurons_extracted = [];  
  
  for j=1:number_of_cells

    % number of spikes by cell
    spike_counts(j) =  fread(fid,1,long_type);
    
    % determine index of specified neurons
    index = find(neurons == cell_ids(j));
    
    % test if in list of specified neurons
    if ~isempty(index)
      
      % debugging/verbose info (uncomment when needed)
      %disp(['loading ' num2str(spike_counts(j)) ' spikes for cell ID ' ...
      %	    num2str(cell_ids(j)) ' (' num2str(j) ')']);

      % actual spike times
      % Note: concatenation is to append multiple records together
      spikes{index} = [spikes{index}; fread(fid, spike_counts(j), spike_time_type)]; 
      
      % append on neuron
      neurons_extracted = [neurons_extracted cell_ids(j)];
      
    else
      % skip over data but do not save it
      fseek(fid, long_size * spike_counts(j), 0);
      %fread(fid, spike_counts(j), spike_time_type);
    end
    
  end  
end

% close file
fclose(fid);

% check that extracted correct number of neurons
if length(neurons_extracted) ~= length(neurons)
  error(['load_rrs_neurons: failed to load all request neurons']);
end

%% 8. Clean Up Output  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract triggers
triggers = spikes{1};

% determine if just triggers
if length(neurons) ~= 1
  % A. Spike times requested

  % find neurons that are not triggers and duplicates
  indices = find(neurons > 0);
  
  % remove trigger and duplicate cells (always)
  neurons = neurons(indices);
  spikes = spikes(indices);
  electrodes = electrodes(indices);
  
  % determine channels associated with neurons
  for i=1:length(neurons)
    index = find(cell_ids == neurons(i));
    electrodes(i) = channels(index);
  end
else
  % B. Spike times not requested

  % find neurons that are not triggers and duplicates
  indices = find(cell_ids > 0);
  
  % remove trigger and duplicate cells (always)
  cell_ids = cell_ids(indices);
  channels = channels(indices);
  spikes = {};
  
  % no spike times requested, so return all cells
  neurons = cell_ids;
  electrodes = channels;
end

% convert spike times to seconds
for i=1:length(spikes)
  spikes{i} = spikes{i} / sampling_frequency;
end

% convert triggers to seconds
triggers = triggers / sampling_frequency;

% compute duration in seconds
duration = nTicks / sampling_frequency;

% create structure
extras = struct('channels',electrodes,'cell_ids',neurons,...
		'triggers',triggers,'duration',duration,...
		'max_contamination', contamination, ...
		'remove_duplicates', remove_duplicates, ...
		'min_rate',min_spikes / duration);

% update user
disp(['extracted ' num2str(length(spikes)) ' cells.']);
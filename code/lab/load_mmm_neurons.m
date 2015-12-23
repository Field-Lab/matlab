function [spikes, extras] = load_mmm_neurons(full_path, neurons)
% LOAD_MMM_NEURONS    loads up old style .neurons file (cells file)
%
% usage:   [spikes] = load_mmm_neurons(neurons_file);
%          [spikes, extras] = load_mmm_neurons(...);
%          [...] = load_mmm_neurons(neurons_file, neurons);
%
% arguments:  neurons-file - full path of .neurons file
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
% were largely based off of EJC's LISP code "read-cell-file" for
% reading in cell files.
%
% examples:  
%     spikes = load_mmm_neurons('/Data/test/data000.neurons');
%     [spikes, extras] = load_mmm_neurons('/Data/test/data000.neurons',[]);
%     spikes = load_mmm_neurons('/Data/test/data000.neurons,[0 2]);
%    
% see also:  LOAD_NEURONS, LOAD_RRS_NEURONS
%
%  2006-07-24 shlens@salk.edu
%  2005-06-21 shlens@salk.edu
%


%% 1. Definitions %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data types
long_type = 'int32';
unsigned_long_type = 'uint32';
float_type = 'float32';
short_type = 'int16';
unsigned_short_type = 'uint16';

% default sampling frequency
sampling_frequency = 20000;  % Hz

% stimulus frame
% Note: this variable is subtracted from the trigger times
% and spike times to properly offset these times with respect to
% the stimulus. The goal of this calculation is to place the first
% trigger time at 1 stimulus frame time prior to 0 sec. A single
% stimulus frame on the latter MMM experiments was 1/120 s.
% BUT, BE CAREFUL. This offset *might* differ for older data sets
% with different recording rigs (e.g. Stanford experiments,
% etc). This code has only been verified for experiments after
% 2002-05-10.
one_stimulus_frame = 1/120; % seconds
one_stimulus_frame = floor(one_stimulus_frame * sampling_frequency);


%% 2. Open File   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 1) && (nargin ~= 2)
  error('load_mmm_neurons: incorrect number of arguments');
end
  
% open file in appropriate format
fid = fopen(full_path,'r','ieee-be');

if (fid == -1)
  error('load_mmm_neurons: file not found');
end


%% 3. Read Header %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file version
file_version = fread(fid,1,long_type);

% not used
% fseek(fid, 4,-1);
% file_info_position = fread(fid,1,long_type);

fseek(fid,12,-1);
record_position = fread(fid,1,long_type);

fseek(fid,16,-1);
channel_position = fread(fid,1,long_type);

% not used
% fseek(fid,64,-1);
% num_spike_files = fread(fid,1,short_type);

fseek(fid,68,-1);
number_of_records = fread(fid,1,short_type);

fseek(fid,70,-1);
number_of_cells = fread(fid,1,short_type);

fseek(fid,record_position,-1);
record_positions = fread(fid,number_of_records,long_type);

% update user
fprintf(1,['Examining ' num2str(number_of_cells) ' cells (MMM v.' ...
      num2str(file_version) ') ... ']);


%% 4. Allocate Memory %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spikes = cell(number_of_cells,1);
triggers = [];


%% 5. Read Records %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:number_of_records
  
  % go to record position
  fseek(fid,record_positions(i),-1);
  
  % not used 
  start_clock = fread(fid,1,unsigned_short_type);
  end_clock = fread(fid,1,unsigned_short_type);

  % number of spikes, triggers
  number_of_triggers = fread(fid,1,long_type);
  total_spikes = fread(fid,1,long_type);
  number_of_spikes = fread(fid,number_of_cells,long_type);
  
  % simple error check
  if (total_spikes ~= sum(number_of_spikes))
    error('load-mmm-cell: spike count is messed up.');
  end
  
  % accumulate triggers
  triggers = [triggers; fread(fid,number_of_triggers,long_type)];
  
  % skip over unused data
  fseek(fid,(4 * number_of_triggers), 0);
  
  % read spike times
  for j=1:number_of_cells
    spikes{j} = [spikes{j}; fread(fid,number_of_spikes(j),long_type)];
  end
  
end


%% 6. Load Channel Info %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab channels information
fseek(fid,channel_position,-1);
channels = fread(fid, number_of_cells, short_type);

% close it up
fclose(fid);


%% 7. Deal with Offset  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please see comment above about one_stimulus_frame

% compute offset (see comments above!)
offset = triggers(1) + one_stimulus_frame;

% properly offset spike and trigger times
triggers = triggers - offset;
for i=1:length(spikes)
  spikes{i} = spikes{i} - offset;
end


%% 8. Clean Up Output % %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert spike times to seconds
for i=1:length(spikes)
  spikes{i} = spikes{i} / sampling_frequency;
end

% convert triggers to seconds
triggers = triggers / sampling_frequency;

% cell ID's are just the order in which they are read out
cell_ids = [0:length(spikes)-1]';

% convert channels to "correct" numbers following
% lisp code "channel-number-to-lead"
corrected_channels = channel_number_to_lead(channels)';

% guess the duration as the maximum spike time
duration = ceil(max(cell2mat(spikes)));

% extract the relevant neurons
% -- this is a bit of a hack
if ~((nargin == 1) || ischar(neurons))
  if ~isempty(neurons)
    % default behavior
    [cell_ids, spikes, corrected_channels] = ...
	extract_neurons(cell_ids, spikes, corrected_channels, neurons);
  else
    % special case
    spikes = {};
  end
end

% create structure
extras = struct('channels',corrected_channels,'cell_ids', cell_ids, ...
		'triggers',triggers,'duration',duration);

% update user
disp([' extracted ' num2str(length(spikes)) ' cells.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ids_out, spikes_out, channels_out] = extract_neurons(...
    ids_in, spikes_in, channels_in, neurons);
% this is a bit of hack to make this function consistent with
% load_rrs_neurons. I am too lazy to do the "correct" solution
% which is not to read in the spike times in the first place.
% However, because the neurons file is often small and not many
% people use this function, I will instead just prune out the 
% user-selected neurons.

spikes_out = cell(length(neurons),1);
channels_out = zeros(length(neurons),1);
ids_out = zeros(length(neurons),1);
bad_neurons = [];

for i=1:length(neurons)
  % find neuron
  index = find(ids_in == neurons(i));

  % find bad neurons
  if isempty(index) 
    bad_neurons = [bad_neurons neurons(i)]; 
    continue;
  end;
  
  % save it
  spikes_out{i} = spikes_in{index};
  ids_out(i) = ids_in(index);
  channels_out(i) = channels_in(index);
end

% bad neurons
if ~isempty(bad_neurons)
  disp(' ');
  error(['load_mmm_neurons: failed to identify ' num2str(bad_neurons)]);
end
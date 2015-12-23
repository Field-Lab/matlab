function new_neurons_file = save_subset_neurons(neurons_raw_file, cell_ids)
% SAVE_SUBSET_NEURONS   Save a subset of neurons
%
% usage: save_subset_neurons(neurons_file, cell_ids);
%
% arguments:  neurons_file - original neurons file
%                 cell_ids - vector of desired cell IDs to be 
%                            saved in new file
%
% shlens 2009-02-12
%%greschner change to vision

% import java classes
import('edu.ucsc.neurobiology.vision.matlab.*');
import('edu.ucsc.neurobiology.vision.io.*');


% grab sampling frequency
sampling_frequency = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% determine new neurons file name
new_neurons_file = parse_new_neurons_name(neurons_raw_file);

% grab/set old header information
old_neurons = NeuronFile(neurons_raw_file);
header = old_neurons.getHeader;
header.version = NeuronFile.INT_VERSION;

% grab triggers
triggers = old_neurons.getTTLTimes;

% load spike times
[spikes, extras] = load_rrs_neurons(neurons_raw_file, cell_ids);

% convert spike times to int32
for i=1:length(cell_ids)
  spikes{i} = int32(spikes{i} * sampling_frequency);
end

% grab electrodes, channels
electrodes = extras.channels;

% save out the neurons file
save_neurons(new_neurons_file, header, triggers, spikes, cell_ids, electrodes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = parse_new_neurons_name(in)

index = findstr(in,'.neurons');
out = [in(1:index-1) '-subset.neurons'];



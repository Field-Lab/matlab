function new_neurons_file = save_split_neurons(neurons_file, cell_ids)
% Split neurons eg shuffel correction STA-correlation bias
%
% arguments:  neurons_file - original neurons file
%                 cell_ids - vector of desired cell IDs to be 
%                            saved in new file
%
%%greschner

% import java classes
import('edu.ucsc.neurobiology.vision.matlab.*');
import('edu.ucsc.neurobiology.vision.io.*');


% grab sampling frequency
sampling_frequency = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% determine new neurons file name
new_neurons_file = parse_new_neurons_name(neurons_file);

% grab/set old header information
old_neurons = NeuronFile(neurons_file);
header = old_neurons.getHeader;
header.version = NeuronFile.INT_VERSION;

% grab triggers
triggers = old_neurons.getTTLTimes;

% load spike times
[spikes, extras] = load_rrs_neurons(neurons_file, cell_ids);

if 0 %rand
    % convert spike times to int32
    for i=1:length(cell_ids)
      sp = int32(spikes{i} * sampling_frequency);
      t=randperm(length(sp));
      spikes_a{i}=sp(sort(t(1:floor(length(t)/2))));
      spikes_b{i}=sp(sort(t(floor(length(t)/2)+1:end)));
    end
end
if 1 %interval
    nr=40;
    gap=1/120*40;   
    for i=1:length(cell_ids)      
        spikes_a{i}=[];
        for ii=1:2:nr  
            t=find(spikes{i}>extras.duration/nr*(ii-1)+gap & spikes{i}<extras.duration/nr*ii);
            spikes_a{i}=[spikes_a{i}; spikes{i}(t)];
        end
        spikes_a{i}=int32(spikes_a{i} * sampling_frequency);
        
        spikes_b{i}=[];
        for ii=2:2:nr  
            t=find(spikes{i}>extras.duration/nr*(ii-1)+gap & spikes{i}<extras.duration/nr*ii);
            spikes_b{i}=[spikes_b{i}; spikes{i}(t)];
        end
        spikes_b{i}=int32(spikes_b{i} * sampling_frequency);        
    end
end

% grab electrodes, channels
electrodes = extras.channels;

% save out the neurons file
save_neurons(new_neurons_file(:,1), header, triggers, spikes_a, cell_ids, electrodes);
save_neurons(new_neurons_file(:,2), header, triggers, spikes_b, cell_ids, electrodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = parse_new_neurons_name(in)

index = findstr(in,'.neurons');
out(:,1) = [in(1:index-1) '-split-1.neurons'];
out(:,2) = [in(1:index-1) '-split-2.neurons'];















function [spikes_weighted,electrodes_used,spike_waveform_points_to_use] = spike_weight_calculator(...
    spikes_original,how_specified,weights,num_electrodes,points_per_electrode);

% figure out which samples to use for PCA

samples_per_spike = size(spikes_original,2)-1;

switch how_specified
    case 'electrodes'
        %interpret "weights" as the list of electrodes to use
        electrodes_used = weights;
        
        % create matrix where each column is an electrode, rows indicate which waveform index is that electrode's data
        temp=reshape(2:num_electrodes*points_per_electrode+1,points_per_electrode,num_electrodes);
        % select subset of points which correspond to the desired electrodes
        spike_waveform_points_to_use=reshape(temp(:,weights),1,points_per_electrode*length(weights));
        
        %delete unused electrodes from the spike waveforms
        spikes_weighted=spikes_original(:,[1 spike_waveform_points_to_use]);
        
    case 'samples'
        %zero any negative values in the weight function
        weight_function = weights.*(weights > 0);
        %multiply each spike by the weight function
        spikes_weighted = spikes_original(:,2:samples_per_spike+1).*(weight_function*ones(samples_per_spike,1));
        spikes_weighted = [spikes_original(:,1) spikes_weighted];
        
        %compute which electrodes are used
        electrodes_used = [];
        for ee = 1:num_electrodes
            min_ = (ee-1)*points_per_electrode+1;
            max_ = min_ + points_per_electrode - 1;
            if max(weights(min_:max_)) > 0  %if it has some > 0 values, include it
                electrodes_used = [electrodes_used ee];
            end
        end
        
        spike_waveform_points_to_use = find(weight_function > 0) + 1;
        
        %return only those electrodes which have at least one used sample
        %[spikes_weighted,electrodes_used,spike_waveform_points_to_use] = spike_weight_calculator(...
        %    spikes_weighted,'electrodes',electrodes_used,num_electrodes,points_per_electrode);
                

end


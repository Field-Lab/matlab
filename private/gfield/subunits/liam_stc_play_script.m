% Get the information for a single cell
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_10.00');
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

% get cones associated with one of the parasol cells
cell_type = 2;
datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);


% get weights for cell type of interest
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.10, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'U');   


plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 30, 'figure', 10)

% get index to RGC of interest
temp_cell_id = 3905;
temp_index = find(datarun.cell_types{cell_type}.cell_ids == temp_cell_id);
temp_cell_index = get_cell_indices(datarun, temp_cell_id);

% get list of cones connected to RGC
temp_cone_ids = find(selection(:,temp_index) == 1);

% check that you have the right set of cones by reconstructing the RF
cone_rfs = Wc(:, temp_cone_ids);
summed_images = sum(cone_rfs, 2);
summed_images = full(summed_images);
reshaped_image = reshape(summed_images, [320,320,3]);
figure(1)
imagesc(reshaped_image)   % plot cones in RF 

% get the cone stimuli for the cell of interest
temp_cone_stim = cone_inputs(:, temp_cone_ids);

% get the spike times for the cell of interest
temp_spike_times = spike_times{temp_cell_index};

% get frames associated with spikes;
spike_frames = temp_cone_stim(temp_spike_times,:);

% make sure cone stimulus is zero mean
cone_stim = temp_cone_stim - repmat(mean(temp_cone_stim, 2), 1,size(temp_cone_stim,2));  
cone_prior = cov(cone_stim);
cone_resp = spike_frames - repmat(mean(spike_frames, 2), 1, size(spike_frames,2));
cone_post = cov(cone_resp);

% Possible computation 1
[eigvecs, eigvals] = eig(cone_prior - cone_post);

% Possible computation 2
[eigvecs, eigvals] = eig(cone_prior^(-0.5) * cone_post * cone_prior^(-0.5) - eye(size(cone_post)));



% project cone RFs through eigenvector 
figure(1)
first_cone_egvec = cone_rfs * eigvecs(:,1);
first_cone_egvec = first_cone_egvec-min(first_cone_egvec);
first_cone_egvec = first_cone_egvec ./ max(first_cone_egvec);
first_cone_egvec = reshape(first_cone_egvec, [320, 320, 3]);
imagesc(first_cone_egvec)
title('PC 1')

% project cone RFs through eigenvector 
figure(2)
first_cone_egvec = cone_rfs * eigvecs(:,2);
first_cone_egvec = first_cone_egvec-min(first_cone_egvec);
first_cone_egvec = first_cone_egvec ./ max(first_cone_egvec);
first_cone_egvec = reshape(first_cone_egvec, [320, 320, 3]);
imagesc(first_cone_egvec)
title('PC 2')

% project cone RFs through eigenvector 
figure(3)
first_cone_egvec = cone_rfs * eigvecs(:,3);
first_cone_egvec = first_cone_egvec-min(first_cone_egvec);
first_cone_egvec = first_cone_egvec ./ max(first_cone_egvec);
first_cone_egvec = reshape(first_cone_egvec, [320, 320, 3]);
imagesc(first_cone_egvec)
title('PC 3')

% project cone RFs through eigenvector 
figure(4)
first_cone_egvec = cone_rfs * eigvecs(:,4);
first_cone_egvec = first_cone_egvec-min(first_cone_egvec);
first_cone_egvec = first_cone_egvec ./ max(first_cone_egvec);
first_cone_egvec = reshape(first_cone_egvec, [320, 320, 3]);
imagesc(first_cone_egvec)
title('PC 4')

% note start and end times, and set time offset
switch datarun.names.nickname
    case 'blueberry'
        start_time = 6400;
        end_time = 9550;
    case 'peach'
        start_time = 0;
        end_time = 2399;
    case 'kiwi'
        start_time = 0;
        end_time = 7440;
    otherwise
        error('start time and end time not set')
end

% generate noise eigenvalues by creating shifting spike times relative to stimulus
num_iters = 100;
time_offsets = -6-round(10000*rand(num_iters,1));
shifted_egvals = zeros(num_iters, length(temp_cone_ids));
for iter = 1:num_iters

    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{temp_cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    shifted_spike_rate_ = circshift(spike_rate_(start_stim:end_stim),time_offsets(iter));

    % translate to spike times (with duplicates for multiple spikes per time bin)
    shifted_spike_times = [];
    for nn = 1:max(spike_rate_)
        shifted_spike_times = [shifted_spike_times; find( shifted_spike_rate_ > (nn-1) )];
    end
    shifted_spike_times = sort(shifted_spike_times);

    shifted_spike_frames = temp_cone_stim(shifted_spike_times,:);

    % make sure cone stimulus is zero mean
    shiftedcone_stim = temp_cone_stim - repmat(mean(temp_cone_stim, 2), 1,size(temp_cone_stim,2));  
    cone_prior = cov(cone_stim);
    cone_resp = spike_frames - repmat(mean(spike_frames, 2), 1, size(spike_frames,2));
    cone_post = cov(cone_resp);

% Possible computation 1
[eigvecs, eigvals] = eig(cone_prior - cone_post);

% Possible computation 2
[eigvecs, eigvals] = eig(cone_prior^(-0.5) * cone_post * cone_prior^(-0.5) - eye(size(cone_post)));


    shifted_egvals(iter,:) = shifted_egvals_'./sum(shifted_egvals_);

end



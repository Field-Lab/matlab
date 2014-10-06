
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, '2008-08-27-5_data003_data003_data003-bayes-msf_70.00');
%datarun = import_single_cone_data(datarun, '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_15.00');
%datarun = import_single_cone_data(datarun, '2008-08-27-5_data003_data003_data003-bayes-msf_70.00');
%datarun = import_single_cone_data(datarun, '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00');

datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);



cell_type = 4;
datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 30, 'figure', 8)


sig_threshold = 0.1;
plot_cell_sampling(datarun, {cell_type}, 'fig_or_axes', 21, 'label', true, 'cone_size', 6, 'thresh', sig_threshold, 'polarity', 1)


[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                        'thresh', sig_threshold, 'radius', [0 inf],...
                                        'polarity', 1, 'contiguity', true, 'scale', 3.0);   


cell_id = 2822; % a nice off midget from peach
cell_id_two = 4067;
%plot_cell_sampling(datarun, cell_id, 'fig_or_axes', 9, 'label', true, 'cone_size', 6, 'thresh', sig_threshold, 'polarity', 1)


cell_pointer = find(datarun.cell_types{cell_type}.cell_ids == cell_id);
cell_pointer_two = find(datarun.cell_types{cell_type}.cell_ids == cell_id_two);
cell_index = get_cell_indices(datarun, cell_id);

% get list of cones connected to RGC
cone_ids = find(selection(:,cell_pointer) == 1);

M_cone_indices = find(datarun.cones.types(cone_ids) == 'M');
L_cone_indices = find(datarun.cones.types(cone_ids) == 'L');
S_cone_indices = find(datarun.cones.types(cone_ids) == 'S');

% check that you have the right set of cones by reconstructing the RF
cone_rfs = Wc(:, cone_ids);
summed_images = sum(cone_rfs, 2);
summed_images = full(summed_images);
reshaped_image = reshape(summed_images, [field_height,field_width,3]);
figure(1)
imagesc(reshaped_image)   % plot cones in RF 

% get the cone stimuli for the cell of interest
cone_stim = cone_inputs(:, cone_ids);

% get the spike times for the cell of interest
cell_spike_times = spike_times{cell_index};

% get frames associated with spikes;
spike_frames = cone_stim(cell_spike_times,:);

% compute the spike triggered variance at each frame:
spike_triggered_var = var(spike_frames, 1,1) - var(cone_stim,1,1);
spike_triggered_mean = mean(spike_frames,1);

figure(25)
clf; hold on
plot(abs(spike_triggered_mean(M_cone_indices)), spike_triggered_var(M_cone_indices), 'go')
plot(abs(spike_triggered_mean(L_cone_indices)), spike_triggered_var(L_cone_indices), 'ro')
plot(abs(spike_triggered_mean(S_cone_indices)), spike_triggered_var(S_cone_indices), 'bo')
hold off

figure(26)
plot(abs(spike_triggered_mean), sqrt(spike_triggered_var), 'ko')

figure(27)
plot(abs(spike_triggered_mean), abs(spike_triggered_mean) ./ sqrt(spike_triggered_var), 'ko')


STA = cone_rfs * spike_triggered_mean';
STE = cone_rfs * spike_triggered_var';


figure(10)
imagesc(norm_image(reshape(STA, [320, 320, 3])))
figure(11)
imagesc(norm_image(reshape(STE, [320, 320, 3])))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate the locations of cones on the grid.

cell_indices = get_cell_indices(datarun, {cell_type});
num_cells = length(cell_indices);

for cll = 1:num_cells
    cell_id = datarun.cell_types{cell_type}.cell_ids(cll); % a nice off midget from peach

    cell_pointer = find(datarun.cell_types{cell_type}.cell_ids == cell_id);

    % get list of cones connected to RGC
    cone_ids = find(selection(:,cell_pointer) == 1);

    % get the cone stimuli for the cell of interest
    cone_stim = cone_inputs(:, cone_ids);

    % get the spike times for the cell of interest
    cell_spike_times = spike_times{cell_indices(cll)};

    % get frames associated with spikes;
    spike_frames = cone_stim(cell_spike_times,:);

    % compute the spike triggered variance and mean
    spike_triggered_var = var(spike_frames, 1,1) - var(cone_stim,1,1);
    spike_triggered_mean = mean(spike_frames,1);

    spike_triggered_var = spike_triggered_var ./ sum(spike_triggered_var);
    spike_triggered_mean = spike_triggered_mean ./ sum(spike_triggered_var);

    figure(1)
    plot(abs(spike_triggered_mean), spike_triggered_var, 'ko')
    hold on
    
    figure(2)
    plot(abs(spike_triggered_mean), abs(spike_triggered_mean) ./ sqrt(abs(spike_triggered_var)), 'ko')
    hold on

    figure(3)
    plot(abs(spike_triggered_mean), sqrt(abs(spike_triggered_var)), 'ko')
    hold on

end














datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_10.00');
%datarun = import_single_cone_data(datarun, '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_15.00');
%datarun = import_single_cone_data(datarun, '2008-08-27-5_data003_data003_data003-bayes-msf_70.00');
%datarun = import_single_cone_data(datarun, '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00');

datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);



cell_type = 2;
datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 30, 'figure', 8)


sig_threshold = 0.05;
plot_cell_sampling(datarun, {cell_type}, 'fig_or_axes', 21, 'label', true, 'cone_size', 6, 'thresh', sig_threshold, 'polarity', 1)


[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                        'thresh', sig_threshold, 'radius', [0 inf],...
                                        'polarity', 1, 'contiguity', true, 'scale', 3.0);   


cell_id_A = 1008;
cell_id_B = 2794;

cell_pointer_A = find(datarun.cell_types{cell_type}.cell_ids == cell_id_A);
cell_pointer_B = find(datarun.cell_types{cell_type}.cell_ids == cell_id_B);
cell_index_A = get_cell_indices(datarun, cell_id_A);
cell_index_B = get_cell_indices(datarun, cell_id_B);

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% cell A
cone_ids_A = find(selection(:,cell_pointer_A) == 1);

% get the cone stimuli for the cell of interest
cone_stim = cone_inputs(:, cone_ids_A);

% get the spike times for the cell of interest
cell_spike_times_A = spike_times{cell_index_A};

% get frames associated with spikes;
spike_frames_A = cone_stim(cell_spike_times_A,:);

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% cell B
cone_ids_B = find(selection(:,cell_pointer_B) == 1);

% get the cone stimuli for the cell of interest
cone_stim = cone_inputs(:, cone_ids_B);

% get the spike times for the cell of interest
cell_spike_times_B = spike_times{cell_index_B};

% get frames associated with spikes;
spike_frames_B = cone_stim(cell_spike_times_B,:);

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% joint spikes

joint_spike_frames = intersect(cell_spike_times_A, cell_spike_times_B);

union_cones = union(cone_ids_A, cone_ids_B);

union_cone_stim = cone_inputs(:,union_cones);

joint_frames = union_cone_stim(joint_spike_frames,:);

[PCs, weights, EgVals] = princomp(joint_frames);

[IDX, VQs] = kmeans(joint_frames,3);
PCs = VQs';

cone_rfs = Wc(:, union_cones);
num_PCs = 3;
for comp = 1:num_PCs;
    figure(comp)
    cone_egvec = cone_rfs * PCs(:,comp);
    cone_egvec = reshape(cone_egvec, [field_height,field_width,3]);
    imagesc(norm_image(cone_egvec))
    title_text = ['PC ', num2str(comp)];
title(title_text)
end

STA_joint = mean(joint_frames,1);
STA_joint = cone_rfs * STA_joint';
figure
imagesc(norm_image(reshape(STA_joint,[320,320,3])))








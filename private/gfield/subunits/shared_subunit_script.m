cell_type = 2;
datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 30, 'figure', 10)


sig_threshold = 0.05;
plot_cell_sampling(datarun, {2}, 'fig_or_axes', 11, 'label', true, 'cone_size', 6, 'thresh', sig_threshold, 'polarity', 1)


[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                        'thresh', sig_threshold, 'radius', [0 inf],...
                                        'polarity', 1, 'contiguity', true, 'scale', 3.0);   




clear reshaped_image_B reshapsed_image_A


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL A
cell_id_A = 1008; % a nice off parasol from peach
cell_pointer_A = find(datarun.cell_types{cell_type}.cell_ids == cell_id_A);
cell_index_A = get_cell_indices(datarun, cell_id_A);

% get list of cones connected to RGC
cone_ids = find(selection(:,cell_pointer_A) == 1);

% check that you have the right set of cones by reconstructing the RF
cone_rfs_A = Wc(:, cone_ids);
summed_images = sum(cone_rfs_A, 2);
summed_images = full(summed_images);
reshaped_image_A = reshape(summed_images, [field_height,field_width,3]);

% get the cone stimuli for the cell of interest
cone_stim = cone_inputs(:, cone_ids);

% get the spike times for the cell of interest
cell_spike_times = spike_times{cell_index_A};

% get frames associated with spikes;
spike_frames = cone_stim(cell_spike_times,:);

% spike_frames is the spike triggered cone-stimulus ensemble
% the next step is to compute the spike triggered average as a check.

[PCs_A, weights_A, eigvals_A] = princomp(spike_frames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL B
cell_id_B = 2794;
cell_pointer_B = find(datarun.cell_types{cell_type}.cell_ids == cell_id_B);
cell_index_B = get_cell_indices(datarun, cell_id_B);

% get list of cones connected to RGC
cone_ids = find(selection(:,cell_pointer_B) == 1);

% check that you have the right set of cones by reconstructing the RF
cone_rfs_B = Wc(:, cone_ids);
summed_images = sum(cone_rfs_B, 2);
summed_images = full(summed_images);
reshaped_image_B = reshape(summed_images, [field_height,field_width,3]);

% get the cone stimuli for the cell of interest
cone_stim = cone_inputs(:, cone_ids);

% get the spike times for the cell of interest
cell_spike_times = spike_times{cell_index_B};

% get frames associated with spikes;
spike_frames = cone_stim(cell_spike_times,:);

% spike_frames is the spike triggered cone-stimulus ensemble
% the next step is to compute the spike triggered average as a check.

[PCs_B, weights_B, eigvals_B] = princomp(spike_frames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualize the two cells and the degree of overlap.
figure(1)
% change the color of the cones feeding B for display purposes

sig_indices = find(squeeze(reshaped_image_B(:,:,1)) > 0);
red_shift = zeros(size(squeeze(reshaped_image_B(:,:,1))));
red_shift(sig_indices) = 0.3;
reshaped_image_B(:,:,1) = reshaped_image_B(:,:,1) + red_shift;


sig_indices = find(squeeze(reshaped_image_A(:,:,3)) > 0);
blue_shift = zeros(size(squeeze(reshaped_image_A(:,:,3))));
blue_shift(sig_indices) = 0.3;
reshaped_image_A(:,:,3) = reshaped_image_A(:,:,3) + blue_shift;
imagesc(norm_image(reshaped_image_A + reshaped_image_B))   % plot cones in RF 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% align eigenvectors to see if subunits are shared.

ymin = 125; ymax = 210; xmin = 140; xmax = 220;
num_PCs = 4;

red_cone_rfs_A = cone_rfs_A;
for cn = 1:size(cone_rfs_A,2);
    tmp_cone_rf = reshape(full(cone_rfs_A(:,cn)), [320,320,3]);
    tmp_cone_rf(:,:,2:3) = 0;
    tmp_cone_rf = reshape(tmp_cone_rf, [],1);
    red_cone_rfs_A(:,cn) = sparse(tmp_cone_rf);
    red_cone_rfs_A(:,cn) = red_cone_rfs_A(:,cn) ./ max(red_cone_rfs_A(:,cn));
end

green_cone_rfs_B = cone_rfs_B;
for cn = 1:size(cone_rfs_B,2);
    tmp_cone_rf = reshape(full(cone_rfs_B(:,cn)), [320,320,3]);
    tmp_cone_rf(:,:,1) = 0;
    tmp_cone_rf(:,:,3) = 0;
    tmp_cone_rf = reshape(tmp_cone_rf, [],1);
    green_cone_rfs_B(:,cn) = sparse(tmp_cone_rf);
    green_cone_rfs_B(:,cn) = green_cone_rfs_B(:,cn) ./ max(green_cone_rfs_B(:,cn));
end

A_PC = 1;
B_PC = 1;

figure(10)
imagesc(norm_image(reshape(red_cone_rfs_A * PCs_A(:,A_PC),[320,320,3]))) 
axis([xmin, xmax, ymin, ymax])
figure(11)
imagesc(norm_image(reshape(green_cone_rfs_B * PCs_B(:,B_PC),[320,320,3]))) 
axis([xmin, xmax, ymin, ymax])

figure(13)
tmp_PC_A = reshape(red_cone_rfs_A * PCs_A(:,A_PC),[320,320,3]);
tmp_PC_B = reshape(green_cone_rfs_B * PCs_B(:,B_PC),[320,320,3]);
imagesc(norm_image(tmp_PC_A + tmp_PC_B));
axis([xmin, xmax, ymin, ymax])

for cmp = 1:num_PCs
    figure(10+cmp)
    imagesc(norm_image(reshape(cone_rfs_A * PCs_A(:,cmp), [320,320,3])))
    axis([xmin, xmax, ymin, ymax])

    figure(20+cmp)
    imagesc(norm_image(reshape(cone_rfs_B * PCs_B(:,cmp), [320,320,3])))
    axis([xmin, xmax, ymin, ymax])

    figure

end












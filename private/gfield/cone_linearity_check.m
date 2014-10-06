% Get the information for a single cell
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, 'peach');
%datarun = import_single_cone_data(datarun, '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_15.00');
%datarun = import_single_cone_data(datarun, '2008-08-27-5_data003_data003_data003-bayes-msf_70.00');
%datarun = import_single_cone_data(datarun, '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00');

datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

%%
% get index to RGC of interest
cell_type = 2;
cell_id = 3334; % a nice off parasol from peach
%cell_id = 3287; % a nice on parasol from peach
%cell_id = 4277;  % a nice off midget from peach
%cell_id = 7051;


% get cone information associated with a cell type
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.1, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0);   


cell_pointer = find(datarun.cell_types{cell_type}.cell_ids == cell_id);
cell_index = get_cell_indices(datarun, cell_id);

% get list of cones connected to RGC
cone_ids = find(selection(:,cell_pointer) == 1);

% check that you have the right set of cones by reconstructing the RF
cone_rfs = Wc(:, cone_ids);
summed_images = sum(cone_rfs, 2);
summed_images = full(summed_images);
baseline= 0.5*ones(field_height, field_width,3);
summed_images = summed_images ./ max(summed_images) ./ 2;
reshaped_image = reshape(summed_images, [field_height,field_width,3]) + baseline;
figure(1)
imagesc(reshaped_image)   % plot cones in RF 
axis([ 110 180 150 220])
axis square
axis off


%%
% neighbor cone indices
cone_a = 63;
cone_b = 61;
cone_pair = [cone_ids(cone_a), cone_ids(cone_b)];

cone_rfs = Wc(:, cone_pair);
summed_images = sum(cone_rfs, 2);
summed_images = full(summed_images);
baseline= 0.5*ones(field_height, field_width,3);
summed_images = summed_images ./ max(summed_images) ./ 2;
reshaped_image = reshape(summed_images, [field_height,field_width,3]) + baseline;
figure
imagesc(reshaped_image)   % plot cones in RF 
axis([ 110 180 150 220])
axis square
axis off

cone_a_inputs = cone_inputs(:,cone_ids(cone_a));
cone_b_inputs = cone_inputs(:,cone_ids(cone_b));
cell_spike_times = spike_times{cell_index};

% expland spike times to correspond to generator signals
num_frames = size(cone_inputs,1);
cell_output = zeros(num_frames,1);
for frm = 1:num_frames
    num_spikes = length(find(cell_spike_times == frm));
    if isempty(num_spikes)
        cell_output(frm) = 0;
    else
        cell_output(frm) = num_spikes;
    end
end

figure(1)
plot(cone_a_inputs, cell_output, 'k.');
figure(2)
plot(cone_b_inputs, cell_output, 'k.');

fit_params_a = fit_static_NL(cell_spike_times, cone_a_inputs)
fit_params_b = fit_static_NL(cell_spike_times, cone_b_inputs)


%%
% distal cone indices
cone_a = 60;
cone_b = 61;
cone_pair = [cone_ids(cone_a), cone_ids(cone_b)];

cone_rfs = Wc(:, cone_pair);
summed_images = sum(cone_rfs, 2);
summed_images = full(summed_images);
baseline= 0.5*ones(field_height, field_width,3);
summed_images = summed_images ./ max(summed_images) ./ 2;
reshaped_image = reshape(summed_images, [field_height,field_width,3]) + baseline;
figure
imagesc(reshaped_image)   % plot cones in RF 
axis([ 110 180 150 220])
axis square
axis off


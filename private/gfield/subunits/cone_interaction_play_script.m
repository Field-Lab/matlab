% script to test for interaction between subunits

% Define clusters of cones "A" and "B." Find frames where the generator signal is near zero 
% (such that there are plenty of frames).  
% Next, defined new groups of cones "C" and "D" which are randomly chosen. Look for a difference
% in the mean spike rate between "A" and "B" clusters and "C" and "D" clusters.

cell_id = 3905;
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id, 'thresh', 0.10, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'U');   

% get indices to cones feeding RGC
cone_ids = find(selection == 1);

% Use the locations and eigenvector structure to define clusters "A" and "B"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define groups from the data

locations = datarun.cones.centers(cone_ids,:);

% Define Unit A
x_inds = find(locations(:,2) > 215 & locations(:,2) < 232);
y_inds = find(locations(:,1) > 108 & locations(:,1) < 122);
A_inds = intersect(x_inds, y_inds);

% Define Unit B
x_inds = find(locations(:,2) > 205 & locations(:,2) < 222);
y_inds = find(locations(:,1) > 94 & locations(:,1) < 107);
B_inds = intersect(x_inds, y_inds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get groups and spikes from sim;
A_inds = find(IDX == 1);
B_inds = find(IDX == 2);

temp_spike_times = spike_times;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get inputs to A unit
A_inputs = sum(cone_inputs(:,A_inds), 2);
B_inputs = sum(cone_inputs(:,B_inds), 2);

wndw = 6;

A_zero_frames = find(A_inputs > wndw);
B_zero_frames = find(B_inputs < -wndw);
cluster_zero_frames = intersect(A_zero_frames, B_zero_frames);
length(cluster_zero_frames)

cluster_spk_rate = length(find(ismember(temp_spike_times, cluster_zero_frames))) ./ (refresh_time * length(cluster_zero_frames));
num_spikes = length(find(ismember(temp_spike_times, cluster_zero_frames)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define random groups of cones with same numbers

% set number of cones in each group
num_grp_C = length(A_inds);
num_grp_D = length(B_inds);

% choose groups
num_iters = 100;
group_spk_rate = zeros(1,num_iters);
for iter = 1:num_iters
    rand('seed', 11111+iter);
    shuffled_cone_indices = randperm(length(cone_ids));
    C_inds = shuffled_cone_indices(1:num_grp_C);
    D_inds = shuffled_cone_indices(num_grp_C+1:num_grp_C+num_grp_D);

    % get the generator potentials across cones in C and D for all stimulus frames
    C_inputs = sum(cone_inputs(:,C_inds),2);
    D_inputs = sum(cone_inputs(:,D_inds),2);

    C_zero_frames = find(C_inputs > wndw);
    D_zero_frames = find(D_inputs < -wndw);
    group_zero_frames = intersect(C_zero_frames, D_zero_frames);

    group_spk_rate(iter) = length(find(ismember(temp_spike_times, group_zero_frames))) ./ (refresh_time * length(group_zero_frames));
    num_spikes = length(find(ismember(temp_spike_times, group_zero_frames)));
end

figure(1); clf
[grp_hist, hist_bins] = hist(group_spk_rate, 20);
bar(hist_bins, grp_hist, 'k')
hold on
plot([cluster_spk_rate, cluster_spk_rate], [0 15], 'r')





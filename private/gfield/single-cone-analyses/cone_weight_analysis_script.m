data_fruit = 'plantain';
radius_scaler = 1.0;
cell_type = 4;
RGC_min_convergence = 5;
RGC_max_convergence = 24;


datarun = import_single_cone_data([], data_fruit);

cell_type_indices = find(datarun.cell_types == cell_type);
[connectivity, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', [],...
                                            'min_convergence', 0, 'max_convergence', 100,...
                                            'max_radius', 2.0, 'normalize', false);
                                        
[center_connectivity, center_connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
                                            'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
                                            'max_radius', radius_scaler, 'normalize', false);
  
[roi_RGC_indices, IA, IB] = intersect(cell_type_indices(connectivity_extras.original_RGC_indices),...
                    cell_type_indices(center_connectivity_extras.original_RGC_indices));



connectivity = connectivity(:,IB);
clear connectivity_extras
                                        
temp_sample_cones = find(connectivity ~= 0);
sampling = connectivity;
sampling(temp_sample_cones) = 1;

connectivity_indices = compute_opponency_index(connectivity, center_connectivity_extras.cone_types);
connectivity_std = std(connectivity_indices)

sampling_indices = compute_opponency_index(sampling, center_connectivity_extras.cone_types);
sampling_std = std(sampling_indices)
                   
[sampling_hist, hist_bins] = hist(sampling_indices, [-1:0.1:1]);
[connectivity_hist, hist_bins] = hist(connectivity_indices, hist_bins);

figure(1)
subplot(2,1,1)
bar(hist_bins, connectivity_hist)
subplot(2,1,2)
bar(hist_bins, sampling_hist)


% generate a list of cone inputs to RGCs
clear temp_num_L temp_num_M
for RGC = 1:center_connectivity_extras.num_RGCs
    temp_cone_indices = find(sampling(:,RGC) == 1);
    temp_cone_ids = center_connectivity_extras.cone_types(temp_cone_indices);
    temp_num_L(RGC) = length(find(temp_cone_ids == 'L'));
    temp_num_M(RGC) = length(find(temp_cone_ids == 'M'));
end
RGC_cone_counts.L = temp_num_L;
RGC_cone_counts.M = temp_num_M;
RGC_cone_counts.difference = temp_num_L - temp_num_M;

hist(RGC_cone_counts.difference,[-15:1:20])
mean(RGC_cone_counts.difference)

% Get RGCs in L dominated regions
L_dom_indices = find(RGC_cone_counts.difference > 5);
% Get RGCs in M dominated regions
M_dom_indices = find(RGC_cone_counts.difference < -1);

num_L_dom = length(L_dom_indices)
num_M_dom = length(M_dom_indices)

% calculate mean CENTER cone weight for each cell in each list
RGC_indices = L_dom_indices;
temp_connectivity = center_connectivity;
L_indices = find(center_connectivity_extras.cone_types == 'L');
M_indices = find(center_connectivity_extras.cone_types == 'M');
clear mean_cone_weight mean_L_weight mean_M_weight
for RGC = 1:length(RGC_indices)
    % get cone weights
    temp_cone_indices = find(temp_connectivity(:,RGC_indices(RGC)) > 0);
    mean_cone_weight(RGC) = mean(temp_connectivity(temp_cone_indices,RGC_indices(RGC)));
    temp_L_indices = intersect(temp_cone_indices, L_indices);
    temp_M_indices = intersect(temp_cone_indices, M_indices);
    if isempty(temp_L_indices)
        mean_L_weight(RGC) = 0;
    else
        mean_L_weight(RGC) = mean(temp_connectivity(temp_L_indices, RGC_indices(RGC)));
    end
    if isempty(temp_M_indices)
        mean_M_weight(RGC) = 0;
    else
        mean_M_weight(RGC) = mean(temp_connectivity(temp_M_indices, RGC_indices(RGC)));
    end
end

M_to_mean_ratio = mean_M_weight ./ mean_cone_weight;
mean(M_to_mean_ratio)

% calculate mean CENTER cone weight for each cell in each list
RGC_indices = M_dom_indices;
temp_connectivity = center_connectivity;
L_indices = find(center_connectivity_extras.cone_types == 'L');
M_indices = find(center_connectivity_extras.cone_types == 'M');
clear mean_cone_weight mean_L_weight mean_M_weight
for RGC = 1:length(RGC_indices)
    % get cone weights
    temp_cone_indices = find(temp_connectivity(:,RGC_indices(RGC)) > 0);
    mean_cone_weight(RGC) = mean(temp_connectivity(temp_cone_indices,RGC_indices(RGC)));
    temp_L_indices = intersect(temp_cone_indices, L_indices);
    temp_M_indices = intersect(temp_cone_indices, M_indices);
    if isempty(temp_L_indices)
        mean_L_weight(RGC) = 0;
    else
        mean_L_weight(RGC) = mean(temp_connectivity(temp_L_indices, RGC_indices(RGC)));
    end
    if isempty(temp_M_indices)
        mean_M_weight(RGC) = 0;
    else
        mean_M_weight(RGC) = mean(temp_connectivity(temp_M_indices, RGC_indices(RGC)));
    end
end

L_to_mean_ratio = mean_L_weight ./ mean_cone_weight;
mean(L_to_mean_ratio)




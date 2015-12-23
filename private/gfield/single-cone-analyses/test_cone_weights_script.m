% Testing Cone weights

cell_type = 4;
radius_scaler = 2;
neighborhoold_radius = 2.0; % identifies cones withing the neighborhood of the RGC, 
                            % not limited to positive cones or cones connected to the cell
data_fruit = 'plantain';
RGC_min_convergence = 7;
RGC_max_convergence = 27;

simulated_mosaic = false;

datarun = import_single_cone_data([], data_fruit);

% exclude S and U cones from type matrix, location matrix, connectivity
% matrix and roi matrix
S_cone_indices = find(datarun.cone_types == 'S');
U_cone_indices = find(datarun.cone_types == 'U');
exclude_cone_indices = [S_cone_indices; U_cone_indices];
keep_cone_indices = setdiff([1:length(datarun.cone_ids)], exclude_cone_indices);
datarun.cone_centers = datarun.cone_centers(keep_cone_indices,:);
datarun.cone_ids = datarun.cone_ids(keep_cone_indices);
datarun.cone_roi = datarun.cone_roi(keep_cone_indices);
datarun.cone_types = datarun.cone_types(keep_cone_indices);
datarun.cone_weights = datarun.cone_weights(keep_cone_indices,:);
datarun

RGC_indices = find(datarun.cell_types == cell_type);
num_RGCs = length(RGC_indices);
num_cones = length(datarun.cone_centers);

new_connectivity_matrix = zeros(num_cones, num_RGCs);
excluded_cell_counter = 0;
for RGC = 1:length(RGC_indices)
    % get number of center cones and associated weights
    temp_RGC_center = datarun.rgc_COMs(RGC_indices(RGC),:);
    if simulated_mosaic == false
        temp_RGC_CA = datarun.rgc_fit_info(RGC_indices(RGC),3) * radius_scaler; % collecting distance for center cones
        temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
        close_cone_indices = find(temp_distances <= temp_RGC_CA);
        positive_cone_indices = find(datarun.cone_weights(:, RGC_indices(RGC)) > 0);
        temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
        temp_num_center_cones = length(temp_cone_indices);
        
        % identify cone types in local neighborhood
        neighborhood_area = datarun.rgc_fit_info(RGC_indices(RGC),3) * neighborhoold_radius;
        cone_neighborhood_indices{RGC} = find(temp_distances <= neighborhood_area);
        %cone_neighborhood_types{RGC} = datarun.cone_types(cone_neighborhood_indices);
    else
        temp_cone_indices = find(datarun.cone_weights(:,RGC) > 0);
        temp_num_center_cones = length(temp_cone_indices);
    end
    
    clear temp_cone_weights
    if temp_num_center_cones < RGC_min_convergence || temp_num_center_cones > RGC_max_convergence
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
    else
        temp_cone_weights = datarun.cone_weights(temp_cone_indices, RGC_indices(RGC));
    end
    
    % assign to new connectivity matrix and normalize
    new_connectivity_matrix(temp_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);    
end
% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(new_connectivity_matrix, [], 1) > 0);
new_connectivity_matrix = new_connectivity_matrix(:, temp_keeper_indices);

[num_roi_cones, num_roi_RGCs] = size(new_connectivity_matrix);

% calculate the OIs
OI_indices = compute_opponency_index(new_connectivity_matrix, datarun.cone_types);
% get a measure across all cells of label purity
non_zero_cone_indices = find(new_connectivity_matrix > 0);
connection_matrix = new_connectivity_matrix;
connection_matrix(non_zero_cone_indices) = 1;
for RGC = 1:num_roi_RGCs
    temp_num_cones = length(find(connection_matrix(:,RGC) == 1));
    connection_matrix(:,RGC) = connection_matrix(:,RGC) ./ temp_num_cones;
end
label_purity_indices = compute_opponency_index(connection_matrix, datarun.cone_types);

% get the extreme OIs
opponent_fraction = 0.25;
nonopponent_fraction = 0.25;
OI_zscores = (abs(zscore(OI_indices))); % convert OIs to zscores
[sorted_zscores, zscore_indices] = sort(OI_zscores, 'ascend');
opponent_sorted_zscore_indices = num_roi_RGCs - round(opponent_fraction * num_roi_RGCs):1:num_roi_RGCs;
nonopponent_sorted_zscore_indices = 1:1:round(nonopponent_fraction * num_roi_RGCs);
opponent_cell_indices = zscore_indices(opponent_sorted_zscore_indices);
nonopponent_cell_indices = zscore_indices(nonopponent_sorted_zscore_indices);

% get connectivity matrix for opponent cells
opponent_connectivity_matrix = new_connectivity_matrix(:, opponent_cell_indices);
nonopponent_connectivity_matrix = new_connectivity_matrix(:, nonopponent_sorted_zscore_indices);

% what is the relative balance of L to M cones in the neighborhood of these populations
for RGC = 1:length(opponent_cell_indices)
    temp_cone_indices = find(new_connectivity_matrix(:,opponent_cell_indices(RGC)) > 0);
    temp_cone_types = datarun.cone_types(temp_cone_indices);
    temp_num_cones = length(temp_cone_indices);
    temp_L_indices = find(temp_cone_types == 'L');
    temp_M_indices = find(temp_cone_types == 'M');
    temp_num_L = length(temp_L_indices);
    temp_num_M = length(temp_M_indices);
    temp_label_index = (temp_num_L - temp_num_M) ./(temp_num_L + temp_num_M);
    opponent_cell_purity(RGC) = temp_label_index;
    op_L_M_coord(RGC,:) = [temp_num_L ./ temp_num_cones, temp_num_M ./ temp_num_cones];
    
    op_oi_coored(RGC,:) = [sum(new_connectivity_matrix(temp_cone_indices(temp_L_indices),opponent_cell_indices(RGC))),...
                            sum(new_connectivity_matrix(temp_cone_indices(temp_M_indices), opponent_cell_indices(RGC)))];
end

for RGC = 1:length(nonopponent_cell_indices)
    temp_cone_indices = find(new_connectivity_matrix(:,nonopponent_cell_indices(RGC)) > 0);
    temp_cone_types = datarun.cone_types(temp_cone_indices);
    temp_num_cones = length(temp_cone_indices);
    temp_L_indices = find(temp_cone_types == 'L');
    temp_M_indices = find(temp_cone_types == 'M');
    temp_num_L = length(temp_L_indices);
    temp_num_M = length(temp_M_indices);
    temp_label_index = (temp_num_L - temp_num_M) ./(temp_num_L + temp_num_M);
    nonopponent_cell_purity(RGC) = temp_label_index;
    nonop_L_M_coord(RGC,:) = [temp_num_L ./ temp_num_cones, temp_num_M ./ temp_num_cones];

    nonop_oi_coored(RGC,:) = [sum(new_connectivity_matrix(temp_cone_indices(temp_L_indices),nonopponent_cell_indices(RGC))),...
                            sum(new_connectivity_matrix(temp_cone_indices(temp_M_indices), nonopponent_cell_indices(RGC)))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what is the relative balance of L to M cones in the neighborhood of these populations
for RGC = 1:length(opponent_cell_indices)
    temp_IDs = find(cone_neighborhood_types{opponent_cell_indices(RGC)} > 0);
    temp_cone_types = datarun.cone_types(temp_cone_indices);
    temp_num_cones = length(temp_cone_indices);
    temp_L_indices = find(temp_cone_types == 'L');
    temp_M_indices = find(temp_cone_types == 'M');
    temp_num_L = length(temp_L_indices);
    temp_num_M = length(temp_M_indices);
    temp_label_index = (temp_num_L - temp_num_M) ./(temp_num_L + temp_num_M);
    opponent_cell_purity(RGC) = temp_label_index;
    op_L_M_coord(RGC,:) = [temp_num_L ./ temp_num_cones, temp_num_M ./ temp_num_cones];
    
    op_oi_coored(RGC,:) = [sum(new_connectivity_matrix(temp_cone_indices(temp_L_indices),opponent_cell_indices(RGC))),...
                            sum(new_connectivity_matrix(temp_cone_indices(temp_M_indices), opponent_cell_indices(RGC)))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_bins = [-1:0.1:1];
figure(1)
subplot(2,1,1)
hist(opponent_cell_purity, temp_bins)
subplot(2,1,2)
hist(nonopponent_cell_purity, temp_bins)

figure(1)
clf
hold on
plot(op_L_M_coord(:,1), op_L_M_coord(:,2), 'r.')
plot(op_oi_coored(:,1), op_oi_coored(:,2), 'k.')

figure(2)
clf
hold on
plot(nonop_oi_coored(:,1), nonop_oi_coored(:,2), 'r.')
%plot([-5:15],zeros(1,21), 'k')
%plot(zeros(1,21),[-5:15], 'k')
%axis([-5 15 -5 15])
plot(nonop_L_M_coord(:,1), nonop_L_M_coord(:,2), 'k.')

purity_category = [zeros(length(opponent_cell_purity), 1), opponent_cell_purity'];
opponency_category = [ones(length(OI_indices(opponent_cell_indices)), 1), OI_indices(opponent_cell_indices)'];

figure(3)
clf
hold on
subplot(1,2,1)
hold on
plot(0, opponent_cell_purity, 'ko')
plot(1, OI_indices(opponent_cell_indices), 'ko')
for RGC = 1:length(opponent_cell_indices)
    plot([0 1], [opponent_cell_purity(RGC), OI_indices(opponent_cell_indices(RGC))], 'r-')
end
axis([-1 2 -1 1])
title('opponent cells')
ylabel('purity index')
xlabel('labels : weights')
subplot(1,2,2)
hold on
plot(0, nonopponent_cell_purity, 'ko')
plot(1, OI_indices(nonopponent_cell_indices), 'ko')
for RGC = 1:length(nonopponent_cell_indices)
    plot([0 1], [nonopponent_cell_purity(RGC), OI_indices(nonopponent_cell_indices(RGC))], 'r-')
end
axis([-1 2 -1 1])
title('nonopponent cells')
xlabel('labels : weights')


% what is the cone expectation?
clear L_dom_RGCs M_dom_RGCs m_cone_expectation l_cone_expectation
imbalance_factor = 1;
M_dom_counter = 0;
L_dom_counter = 0;
for RGC = 1:num_roi_RGCs
    temp_cone_indices = find(new_connectivity_matrix(:,RGC) > 0);
    temp_cone_weights = new_connectivity_matrix(temp_cone_indices, RGC);
    expected_cone_weight = mean(temp_cone_weights);
    temp_num_cones = length(temp_cone_indices);
    
    % compare cone types
    L_cone_indices = find(datarun.cone_types(temp_cone_indices) == 'L');
    M_cone_indices = find(datarun.cone_types(temp_cone_indices) == 'M');
    L_cone_indices_nh = find(datarun.cone_types(cone_neighborhood_indices{RGC}) == 'L');
    M_cone_indices_nh = find(datarun.cone_types(cone_neighborhood_indices{RGC}) == 'M');
    
    
    if isempty(M_cone_indices_nh)
        L_cone_number = temp_num_cones;
        M_cone_number = 0;
    elseif isempty(L_cone_indices_nh)
        M_cone_number = temp_num_cones;
        L_cone_number = 0;
    else
        L_cone_number = length(L_cone_indices_nh);
        M_cone_number = length(M_cone_indices_nh);
    end
    
    if L_cone_number > (M_cone_number+imbalance_factor)
        if M_cone_number > 0
            L_dom_counter = L_dom_counter+1;
            %m_cone_weights = new_connectivity_matrix(temp_cone_indices(L_cone_indices_nh), RGC)
            m_cone_weights = datarun.cone_weights(M_cone_indices_nh, RGC);
            m_cone_expectation(L_dom_counter) = mean(m_cone_weights) / expected_cone_weight;
            L_dom_RGCs(L_dom_counter) = RGC;
        end
    elseif M_cone_number > (L_cone_number+imbalance_factor);
        if L_cone_number > 0;
            M_dom_counter = M_dom_counter +1;
            %l_cone_weights = new_connectivity_matrix(temp_cone_indices(L_cone_indices_nh), RGC);
            L_cone_weights = datarun.cone_weights(L_cone_indices_nh, RGC); 
            l_cone_expectation(M_dom_counter) = mean(l_cone_weights) / expected_cone_weight;
            M_dom_RGCs(M_dom_counter) = RGC;
        end
    end
end

figure(1)
temp_bins = [0:0.2:3];
hist(m_cone_expectation, temp_bins)    
mean(m_cone_expectation) 
standard_error = std(m_cone_expectation) ./ length(m_cone_expectation)

   
figure(2)
temp_bins = [0:0.2:3];
hist(l_cone_expectation, temp_bins)    
mean(l_cone_expectation) 
standard_error = std(l_cone_expectation) ./ length(l_cone_expectation)



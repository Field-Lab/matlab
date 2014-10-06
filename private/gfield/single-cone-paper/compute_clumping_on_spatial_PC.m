function [PC_indices, CIs] = compute_clumping_on_spatial_PC(datarun, cone_indices, PCs)


pixel_scale = 10;

[num_cone, num_PCs] = size(PCs);
 
cone_locations = datarun.cones.centers(cone_indices,:);
distances = ipdm(cone_locations);

CIs = zeros(num_PCs,1);

for pc = 1:num_PCs

    positive_cone_indices = find(PCs(:,pc) > 0);
    negative_cone_indices = find(PCs(:,pc) < 0);

    dists_of_pos_cones = unique(distances(positive_cone_indices, positive_cone_indices));
    non_zero_indices = find(dists_of_pos_cones > 0);
    dists_of_pos_cones = dists_of_pos_cones(non_zero_indices);

    pos_clump_index = sum(exp(-((dists_of_pos_cones./pixel_scale).^2)./2));


    dists_of_neg_cones = unique(distances(negative_cone_indices, negative_cone_indices));
    non_zero_indices = find(dists_of_neg_cones > 0);
    dists_of_neg_cones = dists_of_neg_cones(non_zero_indices);

    neg_clump_index = sum(exp(-((dists_of_neg_cones./pixel_scale).^2)./2));

    CIs(pc) = pos_clump_index + neg_clump_index;

end

[sorted_CIs, indices_for_sorting] = sort(CIs, 'descend');
 
CIs = sorted_CIs;

PC_indices = indices_for_sorting;


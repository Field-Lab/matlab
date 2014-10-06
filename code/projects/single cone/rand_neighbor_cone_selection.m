function new_indexes = rand_neighbor_cone_selection(indexes, connectivity)
% For simulated annealing


% Use this for non-shared cones only
%valid_cones = find(sum(connectivity') == 1);

valid_cones = find(sum(connectivity,2) > 0);

% TODO: This does slightly funny things if there are shared cones...
valid_cones = setdiff(valid_cones, indexes(:));

new_cone = valid_cones(randi(numel(valid_cones)));

possible_rgcs = find(connectivity(new_cone,:));

if length(possible_rgcs) > 1
    rgc = possible_rgcs(randi(length(possible_rgcs)));
else
    rgc = possible_rgcs;
end

cones_per_rgc = size(indexes, 2);
cone = randi(cones_per_rgc);
new_indexes = indexes;
new_indexes(rgc,cone) = new_cone;
function fd = clump_scale_constant(cones)

locations = [cones(:,2) cones(:,3)];
distances = squareform(pdist(locations));
temp = [];
for d = 1:length(cones(:,1))
    temp = [temp min(distances(d, distances(d,:) > 0))]; %#ok<AGROW>
end

nearest = temp;
fd = median(nearest);


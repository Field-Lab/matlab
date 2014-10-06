function spaced = space_voronoi_masks(masks, centers, near_cones_list, neighbor_min_dist, self_max_dist)
% SPACE_VORONOI_MASKS   Usually called as DATARUN (@STRUCT) version
% usage: spaced = space_voronoi_masks(masks, centers, near_cones_list, neighbor_min_dist, self_max_dist)
%
% 2012-09 phli, broken out of datarun version
%

neighbor_min_dist_squared = neighbor_min_dist * neighbor_min_dist;
self_max_dist_squared     = self_max_dist     * self_max_dist;

spaced = cell(size(masks));
for i = 1:length(spaced)
    mask = masks{i};
    if isempty(mask); continue; end
    
    % Get coordinates for mask pixels
    [r,c] = find(mask);
    
    % Run through mask pixels and unset them if they are too close to other cones
    if neighbor_min_dist > 0
        near_cones = near_cones_list{i};
        near_centers = centers(near_cones,:);        
        for j = 1:length(r)
            delx = near_centers(:,1) - c(j);
            dely = near_centers(:,2) - r(j);
            distsqr = delx.^2 + dely.^2;
            
            if any(distsqr < neighbor_min_dist_squared)
                mask(r(j),c(j)) = false;
            end
        end
    end
    
    % Run through mask pixels and unset them if they are too far from own cone
    if self_max_dist < Inf
        this_center = centers(i,:);
        for j = 1:length(r)
            delx = this_center(1) - c(j);
            dely = this_center(2) - r(j);
            distsqr = delx.^2 + dely.^2;
            
            if distsqr > self_max_dist_squared
                mask(r(j),c(j)) = false;
            end
        end
    end
        
    spaced{i} = mask;
end
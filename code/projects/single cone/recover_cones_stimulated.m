function cones_by_index = recover_cones_stimulated(conerun, map)
% RECOVER_CONES_STIMULATED
% usage: cones_by_index = recover_cones_stimulated(conerun, map)
%
% We weren't good about writing down which cones were used to generate cone
% stim maps.  This will attempt to recover the cone stimulus settings based 
% on the map value over the cone center.
%
% 2012, phli
%

conesx = conerun.cones.centers(:,1) .* conerun.stimulus.stixel_width;
conesy = conerun.cones.centers(:,2) .* conerun.stimulus.stixel_height;
conesi = sub2ind(size(map), round(conesy), round(conesx));
cone_mapvals = map(conesi);

uindices = setdiff(unique(cone_mapvals), 0);
cones_by_index = cell(1, max(map(:)));
for i = 1:length(uindices)
    index = uindices(i);
    cones_by_index{index} = find(cone_mapvals == index);
end
function [full_sta, cones] = expand_voronoi_sta(raw_sta, map)

max_cone = max(map(:));

full_sta=zeros(size(map,1),size(map,2),size(raw_sta,2));
cones = zeros(max_cone,2);

for i=1:max_cone
    [a, b] = find(map==i);
  
    if ~isempty(a)
        cones(i,1) = mean(a);
        cones(i,2) = mean(b);
        for j = 1:length(a)
            full_sta(a(j),b(j),:) = raw_sta(i,:);
        end
    else
        cones(i,:) = nan;
    end
end

function [full_sta, cones] = expand_voronoi_sta(raw_sta, map)

max_cone = max(map(:));

full_sta=zeros(size(map,1),size(map,2),size(raw_sta,2));
cones = zeros(max_cone,2);

for i=1:max_cone-1
    [a, b] = find(map==i);
  
    if ~isempty(a)
        cones(i+1,1) = mean(a);
        cones(i+1,2) = mean(b);
        for j = 1:length(a)
            full_sta(a(j),b(j),:) = raw_sta(i+1,:);
        end
    else
        cones(i+1,:) = nan;
    end
end
% next part is VERY questionable
i = max_cone;
[a, b] = find(map==i);
cones(1,1) = mean(a);
cones(1,2) = mean(b);
for j = 1:length(a)
    full_sta(a(j),b(j),:) = raw_sta(1,:);
end

%     figure
%     imagesc(full_sta)
%     hold on
%     for i=1:569
%         text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
%     end


%   center_cones = find(raw_sta(:,27)<=-0.05)-1;
%     center_cones(isnan(cones(center_cones,1)))=[];
%     figure
%     imagesc(full_sta(:,:,27))
%     hold on
%     for i=center_cones
%         text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
%     end

function plot_maps_for_interactions(center_cones, map, cones, full_sta, datarunID, save_path, sta1, sta2, cones1, cones2, conts)


path2save = [save_path, int2str(datarunID)];
if ~isdir(path2save)
    mkdir(path2save);
    mkdir([path2save, '/tiff/']);
    mkdir([path2save, '/svg/']);
end

bord = 20;
coords = cones(center_cones,:);
y = [min(coords(:,1))-bord max(coords(:,1))+bord];
x = [min(coords(:,2))-bord max(coords(:,2))+bord];
if diff(x)>diff(y)
    d = (diff(x)-diff(y))/2;
    y(1) = y(1)-d;
    y(2) = y(2)+d;
else
    d = (diff(y)-diff(x))/2;
    x(1) = x(1)-d;
    x(2) = x(2)+d;
end

% for ej yale talk, cell 3736 - matching Jeremy coordinates
% figure
% subplot('position', [0 0 1 1])
% colormap gray
% imagesc(full_sta)
% axis([302.5 344.5  476.5 502.5])
% set(gca, 'dataaspectratio', [1 1 1], 'visible', 'off')
% hold on
% for i=1:length(center_cones)
%     plot(conts{center_cones(i),1}, conts{center_cones(i),2}, 'color', [1 0 0], 'linewidth', 0.2)
% end



figure
set(gcf, 'position', [422  140  975  944])

subplot(2,2,1)

tmp =sta1/max(abs(sta1(:)))/2+0.5;
tmp_map = map;
tmp_map(map>0.5) = 0.1;
tmp = repmat(tmp, 1, 1, 3);
tmp(:,:,1) = tmp(:,:,1) + tmp_map;
imagesc(tmp)
hold on
for i=1:length(center_cones)
    plot(conts{center_cones(i),1}, conts{center_cones(i),2}, 'color', [1 0 0], 'linewidth', 0.4)
end
plot(cones1(:,1), cones1(:,2), 'x', 'markersize', 3,'color', [1 1 0])
plot(cones2(:,1), cones2(:,2), '+', 'markersize', 3,'color', [0 1 1])
title('single cone run 1')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

subplot(2,2,2)
tmp =sta2/max(abs(sta2(:)))/2+0.5;
tmp_map = map;
tmp_map(map>0.5) = 0.1;
tmp = repmat(tmp, 1, 1, 3);
tmp(:,:,1) = tmp(:,:,1) + tmp_map;
imagesc(tmp)
hold on
for i=1:length(center_cones)
    plot(conts{center_cones(i),1}, conts{center_cones(i),2}, 'color', [1 0 0], 'linewidth', 0.3)
end
plot(cones1(:,1), cones1(:,2), 'x', 'markersize', 3,'color', [1 1 0])
plot(cones2(:,1), cones2(:,2), '+', 'markersize', 3,'color', [0 1 1])
title('single cone run 2')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

subplot(2,2,3)
tmp_map = map;
tmp_map(map>0.5) = 0.4;
tmp_map(1,1) = 1;
colormap gray
imagesc(tmp_map)
title('single cone run 1')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
hold on
for i=1:length(cones)
    plot(conts{i,1}, conts{i,2}, 'color', [0 0 1], 'linewidth', 0.2)
end
for i=1:length(center_cones)
    plot(conts{center_cones(i),1}, conts{center_cones(i),2}, 'color', [1 0 0], 'linewidth', 0.3)
end
plot(cones1(:,1), cones1(:,2), 'x', 'markersize', 3,'color', [1 1 0])
plot(cones2(:,1), cones2(:,2), '+', 'markersize', 3,'color', [0 1 1])

subplot(2,2,4)
colormap gray
imagesc(full_sta)
title('voronoi')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
hold on
for i=1:length(center_cones)
    plot(conts{center_cones(i),1}, conts{center_cones(i),2}, 'color', [1 0 0], 'linewidth', 0.2)
    text(min(conts{center_cones(i),1}), min(conts{center_cones(i),2}), int2str(center_cones(i)), 'color', 'g', 'fontsize', 6)
end

%  drawnow
% saveas(gcf,[path2save,'/tiff/ID_',int2str(cellID),'_maps.tiff'])
% saveas(gcf,[path2save,'/svg/ID_',int2str(cellID),'_maps.svg'])
drawnow
close(gcf)

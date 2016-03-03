function offline_contrast_response(inputs, spikes, nbins_cone1, nbins_cone2, center_cones, map,...
    cones, comb, datarunID, cellID,  save_path, sta1, sta2, voronoi_contours, pol)

path2save = [save_path, int2str(datarunID)];
if ~isdir(path2save)
    mkdir(path2save);
    mkdir([path2save, '/tiff/']);
    mkdir([path2save, '/svg/']);
end


spikes(spikes>size(inputs,2)) = []; 
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

size_cone1_bin = floor(size(inputs,2)/nbins_cone1);
size_cone2_bin = floor(size(inputs,2)/nbins_cone2);

colors = [1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3;  0.1 0.1 0.1; ...
    1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3; 0.1 0.1 0.1; ...
    1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3;  0.1 0.1 0.1];
colors2 = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1]/2;

coords = cones(center_cones,:);
y = [min(coords(:,1))-5 max(coords(:,1))+5];
x = [min(coords(:,2))-5 max(coords(:,2))+5];
if diff(x)>diff(y)
    d = (diff(x)-diff(y))/2;
    y(1) = y(1)-d;
    y(2) = y(2)+d;
else
    d = (diff(y)-diff(x))/2;
    x(1) = x(1)-d;
    x(2) = x(2)+d;
end

figure
set(gcf, 'position', [422  140  975  944])

subplot(2,2,1)
colormap gray
imagesc(sta1)
hold on
cnt = 1;cnt1 = 1;
for i=1:max(map(:))
    if ~isempty(voronoi_contours{i, 1})
        if ~isempty(find(center_cones==i+1,1))
            wid = 1.2; col = colors(cnt1, :);
            cnt1=cnt1+1;
        else
            wid = 0.1;
            col = colors2(cnt,:);
            if cnt==5;
                cnt = 0;
            end
            cnt = cnt+1;
        end
        plot(voronoi_contours{i, 1}, voronoi_contours{i, 2},'color', col, 'linewidth', wid)
    end
end
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
title('single cone run 1')

subplot(2,2,2)
colormap gray
imagesc(sta2)
hold on
cnt = 1;cnt1 = 1;
for i=1:max(map(:))
    if ~isempty(voronoi_contours{i, 1})
        if ~isempty(find(center_cones==i+1,1))
            wid = 1.2; col = colors(cnt1, :);
            cnt1=cnt1+1;
        else
            wid = 0.1;
            col = colors2(cnt,:);
            if cnt==5;
                cnt = 0;
            end
            cnt = cnt+1;
        end
        plot(voronoi_contours{i, 1}, voronoi_contours{i, 2},'color', col, 'linewidth', wid)
    end
end
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
title('single cone run 2')

subplot(2,2,3)
imagesc(comb)
hold on
for i=center_cones
    text(cones(i,2)-1,cones(i,1), int2str(i), 'color', 'r', 'fontsize', 5)
end
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
title('Voronoi convolved STA')

subplot(2,2,4)
colored_cones = zeros([size(comb), 3]);
for i = 1:length(center_cones)
    [a, b] = find(map==center_cones(i)-1);
    for j = 1:length(a)
        colored_cones(a(j),b(j), :) = colors(i,:);
    end
end
imagesc(colored_cones);
hold on
for i=center_cones
    text(cones(i,2)-1,cones(i,1), int2str(i), 'color', [1 1 1], 'fontsize', 5)
end
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])
title('Voronoi colored regions')

saveas(gcf,[path2save,'/tiff/ID_',int2str(cellID),'_maps.tiff'])
saveas(gcf,[path2save,'/svg/ID_',int2str(cellID),'_maps.svg'])
close(gcf)


    
for cone1 = 1:length(center_cones)
    
    figure(cone1)
    set(gcf, 'position', [178 123 1208 979])
    
    % referent cone unconditional response
    tmp = sort(inputs(cone1,:));    
    cone1_contr = [tmp(1:size_cone1_bin:end) tmp(end)];
    cone1_uncond_rate = zeros(1, nbins_cone1);
    cone1_valid = cell(1,nbins_cone1);
    for j = 1:nbins_cone1
        cone1_valid{j} = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
        cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
    end
    cone1_xrange = cone1_contr(1:end-1)+diff(cone1_contr)/2;
    
    % conditional response on all other cones
    other_cones = 1:size(inputs,1);
    other_cones(cone1) = [];
    
    pos = [0.02 0.68 0.3 0.28; 0.34 0.68 0.3 0.28; 0.7 0.68 0.28 0.28;...
        0.02 0.35 0.3 0.28; 0.34 0.35 0.3 0.28; 0.7 0.35 0.28 0.28;...
        0.02 0.02 0.3 0.28; 0.34 0.02 0.3 0.28; 0.7 0.02 0.28 0.28];
%     subplot('position', pos(1,:))
    for i=1:9
        h(i) = subplot('position', pos(i,:));
        hold on
    end

    for cone2 = other_cones        
        tmp = sort(inputs(cone2,:));
        cone2_contr = [tmp(1:size_cone2_bin:end) tmp(end)];        
        cond_rate = zeros(nbins_cone2, nbins_cone1-1);
        % highest and lowest contrasts on the other cone ONLY
        for i=[1,nbins_cone2]         
            cone2_valid = find(inputs(cone2,:)>=cone2_contr(i) & inputs(cone2,:)<=cone2_contr(i+1));
            for j = 1:nbins_cone1
                cond_rate(i,j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
            end
            [params{cone2,i} params1{cone2,i} params2{cone2,i} resp{cone2,i} resp1{cone2,i} x1{cone2,i}, resp2{cone2,i}, x2{cone2,i}] =...
                fit_cone_interactions(cone1_xrange, cone1_uncond_rate, cond_rate(i,:));
        end
        subplot(h(1))
        plot(cone1_xrange,resp{cone2,1},'color',colors(cone2,:), 'linewidth', 1)
        subplot(h(2))
        plot(cone1_xrange,resp{cone2,nbins_cone2},'color',colors(cone2,:), 'linewidth', 1)
        subplot(h(4))
        plot(x1{cone2,1},resp1{cone2,1},'color',colors(cone2,:), 'linewidth', 1)
        subplot(h(5))
        plot(x1{cone2,nbins_cone2},resp1{cone2,nbins_cone2},'color',colors(cone2,:), 'linewidth', 1)
        subplot(h(7))
        plot(x2{cone2,1},resp2{cone2,1},'color',colors(cone2,:), 'linewidth', 1)
        subplot(h(8))
        plot(x2{cone2,nbins_cone2},resp2{cone2,nbins_cone2},'color',colors(cone2,:), 'linewidth', 1)
    end
    
    ylims = [];xlims = [];
    
    for i=[1 2 4 5 7 8]
        subplot(h(i))
        axis tight;
        ylims = [ylims get(gca, 'YLim')];
        xlims = [xlims get(gca, 'XLim')];
    end
    ylims = [min(ylims) max(ylims)];
    xlims = [min(xlims) max(xlims)];
    
    for i=[1 2 4 5 7 8]
        subplot(h(i))
        plot(cone1_xrange,cone1_uncond_rate,  '-x', 'color',[0.05 0.05 0.05], 'linewidth',1.5, 'markersize', 3)
        axis([xlims(1)*1.05  xlims(end)*1.05 ylims(1)*0.95 ylims(2)*1.05])        
        set(gca, 'xtick', 0, 'xticklabel', '')
        
        if i==2 || i==5 || i == 8
            tit1 = 'max contrast; ';
            tmpy = [round((ylims(1)+diff(ylims)/3)*100)/100 round((ylims(1)+diff(ylims)/3*2)*100)/100];
            tmpy = unique(tmpy);
            set(gca, 'ytick',tmpy , 'fontsize', 6, 'fontweight', 'normal')
        else
            tit1 = 'min contrast; ';
            set(gca, 'ytick', sort([0 round(ylims(2)*100)/100]), 'yticklabel', '')
        end
        if i<3
            tit='y shifted';
        elseif i<6
            tit='x shifted';
        else            
            tit='x AND y shifted';
        end
        title([tit1 tit], 'fontsize', 8, 'fontweight', 'normal')
    end
    
  
    subplot(h(3))
    clear t
    for i=1:size(params,1)
        if ~isempty(params{i,1})
            t(i,1) = params{i,1}(end);
            t(i,2) = params1{i,1}(end);
            t(i,3) = params2{i,1}(end);
        end
    end
    colored_cones = zeros([size(comb), 3]);
    for i = 1:length(center_cones)
        [a, b] = find(map==center_cones(i));
        if i~=cone1
            c = t(i,:)/max(t(i,:));
        else c= [1 1 1];
        end
        for j = 1:length(a)
            colored_cones(a(j),b(j), :) = c;
        end
    end
    imagesc(colored_cones);
    axis([x y])
    set(gca, 'xtick', [],  'xticklabel', '')
    set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')    
    axis ij
    set(gca, 'dataaspectratio', [1 1 1])
    title({'min contrast; red: x shift better','green: y shift better, white: ref cone'}, 'fontsize', 8, 'fontweight', 'normal')

    subplot(h(6))
    clear t
    for i=1:size(params,1)
        if ~isempty(params{i,1})
            t(i,1) = params{i,end}(end);
            t(i,2) = params1{i,end}(end);
            t(i,3) = params2{i,end}(end);
        end
    end
    colored_cones = zeros([size(comb), 3]);
    for i = 1:length(center_cones)
        [a, b] = find(map==center_cones(i));
        if i~=cone1
            c = t(i,:)/max(t(i,:));
        else c= [1 1 1];
        end
        for j = 1:length(a)
            colored_cones(a(j),b(j), :) = c;
        end
    end
    imagesc(colored_cones);
    axis([x y])
    set(gca, 'xtick', [],  'xticklabel', '')
    set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
    axis ij
    set(gca, 'dataaspectratio', [1 1 1])
    title({'max contrast; red: x shift better','green: y shift better, white: ref cone'}, 'fontsize', 8, 'fontweight', 'normal')
    
    subplot(h(9))
    colormap gray
    imagesc(sta1)
    hold on
    cnt = 1;cnt1 = 1;
    for i=1:max(map(:))
        if ~isempty(voronoi_contours{i, 1})
            if ~isempty(find(center_cones==i,1))
                wid = 1.2; col = colors(cnt1, :);
                cnt1=cnt1+1;
            else
                wid = 0.1;
                col = colors2(cnt,:);
                if cnt==5;
                    cnt = 0;
                end
                cnt = cnt+1;
            end
            plot(voronoi_contours{i, 1}, voronoi_contours{i, 2},'color', col, 'linewidth', wid)
        end
    end
    hold on
    if pol>0
        cols = [0 0 0];
    else
        cols = [1 1 1];
    end
    for i=center_cones
        text(cones(i,2)-1,cones(i,1), int2str(i), 'color', cols, 'fontsize', 5)
    end
    axis([x y])
    set(gca, 'xtick', round([(x(1)+diff(x)/3) (x(1)+diff(x)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
    set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')  
    set(gca, 'dataaspectratio', [1 1 1])
    axis ij
    title('single cone run 1',  'fontsize', 8, 'fontweight', 'normal')


    drawnow
    saveas(gcf,[path2save,'/tiff/', 'ID_',int2str(cellID), '_refcone_' int2str(center_cones(cone1)),'.tiff'])
    saveas(gcf,[path2save,'/svg/', 'ID_',int2str(cellID), '_refcone_' int2str(center_cones(cone1)),'.svg'])
    close(gcf)
end



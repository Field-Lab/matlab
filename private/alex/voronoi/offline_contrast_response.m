function offline_contrast_response(inputs, spikes, nbins_cone1, nbins_cone2, center_cones, map,...
    cones, full_sta, datarunID, cellID,  save_path, sta1, sta2, cones1, cones2, voronoi_contours, pol)

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

colors = [1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3;  0.1 0.1 0.1; ...
    1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3; 0.1 0.1 0.1; ...
    1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0; 0 0.6 0; 1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3;  0.1 0.1 0.1];
colors2 = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1]/2;

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

figure
set(gcf, 'position', [422  140  975  944])

subplot(2,2,1)
comb = zeros(size(sta1, 1), size(sta1, 2), 3);
comb(:,:,1) = map;
comb(comb>0.5) = 0.5;
comb(:,:,2) = sta1/max(abs(sta1(:)));
imagesc(comb)
hold on
plot(cones1(:,1), cones1(:,2), 'x', 'color', [1 1 1])
plot(cones2(:,1), cones2(:,2), '+', 'color', [1 1 0])
title('single cone run 1')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

subplot(2,2,2)
comb = zeros(size(sta1, 1), size(sta1, 2), 3);
comb(:,:,1) = map;
comb(comb>0.5) = 0.5;
comb(:,:,3) = sta2/max(abs(sta2(:)));
imagesc(comb)
hold on
plot(cones1(:,1), cones1(:,2), 'x', 'color', [1 1 1])
plot(cones2(:,1), cones2(:,2), '+', 'color', [1 1 0])
title('single cone run 2')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

subplot(2,2,3)
comb = zeros(size(sta1, 1), size(sta1, 2), 3);
comb(:,:,1) = map;
comb(comb>0.5) = 0.5;
comb(:,:,2) = sta1/max(abs(sta1(:)));
comb(:,:,3) = sta2/max(abs(sta2(:)));
imagesc(comb)
hold on
plot(cones1(:,1), cones1(:,2), 'x', 'color', [1 1 1])
plot(cones2(:,1), cones2(:,2), '+', 'color', [1 1 0])
title('single cone combo')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

subplot(2,2,4)
colormap gray
imagesc(full_sta)
title('voronoi')
axis([x y])
set(gca, 'dataaspectratio', [1 1 1])

% saveas(gcf,[path2save,'/tiff/ID_',int2str(cellID),'_maps.tiff'])
% saveas(gcf,[path2save,'/svg/ID_',int2str(cellID),'_maps.svg'])
% close(gcf)

cone1 =1;cone2 = 2;
nbins_cone1 = 8;
nbins_cone2 = 8;

% 
% tmp = sort(inputs(cone1,:));
% clear cone1_contr
% for j=1:nbins_cone1-1
%     cone1_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone1));
% end
% x = [tmp(1); cone1_contr'; tmp(end)];
% x = x(1:end-1)+diff(x)/2;
% 
% cone1_contr = [-10 cone1_contr 10];
% cone1_uncond_rate = zeros(1, nbins_cone1);
% cone1_valid = cell(1,nbins_cone1);
% for j = 1:nbins_cone1
%     cone1_valid{j} = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
%     cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
% end
% x2 = [x x x];
% 
% tmp = sort(inputs(cone2,:));
% cone2_valid = find(inputs(cone2,:)>=tmp(round(size(tmp,2)*(nbins_cone1-1)/nbins_cone1)));
% cond_rate = zeros(nbins_cone1,1);
% for j = 1:nbins_cone1
%     cond_rate(j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
% end
% y2 = [cone1_uncond_rate' cond_rate];
% 
% cone2_valid = find(inputs(cone2,:)<tmp(round(size(tmp,2)*1/nbins_cone1)));
% cond_rate = zeros(nbins_cone1,1);
% for j = 1:nbins_cone1
%     cond_rate(j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
% end
% y2 = [y2 cond_rate];
% 
% [p g resnorm residual resnorm1 residual1] = fit_2_curves(x2,y2);
% 
marks = 'x+o';

figure
plot(x2(:,1), y2(:,1), '-x')
hold on
plot(x2(:,2), y2(:,2), '-+')
plot(x2(:,3), y2(:,3), '-o')

figure
hold on
for i=1:size(x2,2)
    y = p(end-2) .* normcdf(x2(:,i) + p(i), p(end), p(end-1));
    plot(x2(:,i) + p(i), y, marks(i))
end

figure
hold on
for i=1:size(x2,2)
    y = g(end-2) .* normcdf(x2(:,i), g(end), g(end-1))+g(i);
    plot(x2(:,i), y, marks(i))
end


% 
% dat = [cone1_uncond_rate' cond_rate];
for cone1 = 1:length(center_cones)
    
    nbins_cone1 = 20;
    nbins_cone2 = 8;
    % referent cone unconditional response
    tmp = sort(inputs(cone1,:));
    clear cone1_contr
    for j=1:nbins_cone1-1
        cone1_contr(j) = tmp(round(size(tmp,2)*j/nbins_cone1));
    end
    x = [tmp(1); cone1_contr'; tmp(end)];
    x = x(1:end-1)+diff(x)/2;
    
    
    cone1_contr = [-10 cone1_contr 10];
    cone1_uncond_rate = zeros(1, nbins_cone1);
    cone1_valid = cell(1,nbins_cone1);
    for j = 1:nbins_cone1
        cone1_valid{j} = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
        cone1_uncond_rate(j) = mean(spike_rate(cone1_valid{j}));
    end
    
    xt=[0.1 0.1 0.1];
    [p, ress] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),xt, x, cone1_uncond_rate');
    [g, ress] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),p, inputs(cone1,:)', spike_rate);
    
    x2 = [x x];
    y2 = [cone1_uncond_rate' cond_rate];
    [g, ress] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),p, x2, y2);
    
    
    
    
    
    
    
    xt=[0];
    [k, ress] = lsqcurvefit(@(t,x) normcdf(x,g(1),g(2))+g(3)+t(1),xt, inputs(cone1,cone2_valid)', spike_rate(cone2_valid));
    params = [p k ress];
    resp = spike_rate(cone2_valid)-k(1);
    
    
    xt=[0,1];
    [k, ress] = lsqcurvefit(@(t,x) normcdf(x+t(1),g(1),g(2))+g(3),xt, inputs(cone1,cone2_valid)', spike_rate(cone2_valid), [(x(1)-x(end))*2 -Inf]);
    params1 = [p k ress];
    resp1 = spike_rate(cone2_valid);
    x1 = x+k(1);
    
    
    xt=[0,0];
    [k, ress] = lsqcurvefit(@(t,x) normcdf(x+t(1),g(1),g(2))+g(3)+t(2),xt, inputs(cone1,cone2_valid)', spike_rate(cone2_valid), [(x(1)-x(end))*2 -Inf]);
    params2 = [p k ress];
    resp2 = spike_rate(cone2_valid)-k(2);
    x2 = x+k(1);
    
    
    figure
    subplot(1,3,1)
    plot(inputs(cone1,cone2_valid),resp,'color',colors(cone2,:), 'linewidth', 1)
    subplot(1,3,2)
    plot(x1,resp1,'color',colors(cone2,:), 'linewidth', 1)
    subplot(1,3,3)
    plot(x2,resp2,'color',colors(cone2,:), 'linewidth', 1)
    
    
    
    
    cone2 = 1;
    tmp = sort(inputs(cone2,:));
    cone2_valid = find(inputs(cone2,:)>=tmp(round(size(tmp,2)*(nbins_cone1-1)/nbins_cone1)));
    [q] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),g, inputs(cone1,cone2_valid)', spike_rate(cone2_valid));
    td = inputs(cone1,cone2_valid)'+inputs(cone2,cone2_valid)';
    [k] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),g, td, spike_rate(cone2_valid));
    
    
    % conditional response on all other cones
    other_cones = 1:size(inputs,1);
    other_cones(cone1) = [];
    
    figure
    for cone2 = other_cones
        tmp = sort(inputs(cone2,:));
        cone2_valid = find(inputs(cone2,:)>=tmp(round(size(tmp,2)*(nbins_cone1-1)/nbins_cone1)));
        cond_rate = zeros(nbins_cone1,1);
        for j = 1:nbins_cone1
            cond_rate(j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
        end
        plot(cond_rate)
        hold on
    end
    plot(cone1_uncond_rate, 'color', 'k', 'linewidth', 2)
    
    [params{cone2,i} params1{cone2,i} params2{cone2,i} resp{cone2,i} resp1{cone2,i} x1{cone2,i}, resp2{cone2,i}, x2{cone2,i}] =...
        fit_cone_interactions(cone1_xrange, cone1_uncond_rate, cond_rate);
    
    
    % highest and lowest contrasts on the other cone ONLY
%     for i=[1,nbins_cone2]
%         cone2_valid = find(inputs(cone2,:)>=cone2_contr(i) & inputs(cone2,:)<=cone2_contr(i+1));
%         for j = 1:nbins_cone1
%             cond_rate(i,j) = mean(spike_rate(intersect(cone1_valid{j}, cone2_valid)));
%         end
%         [params{cone2,i} params1{cone2,i} params2{cone2,i} resp{cone2,i} resp1{cone2,i} x1{cone2,i}, resp2{cone2,i}, x2{cone2,i}] =...
%             fit_cone_interactions(cone1_xrange, cone1_uncond_rate, cond_rate(i,:));
%     end
%     
    
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
% 
% ylims = [];xlims = [];
% 
% for i=[1 2 4 5 7 8]
%     subplot(h(i))
%     axis tight;
%     ylims = [ylims get(gca, 'YLim')];
%     xlims = [xlims get(gca, 'XLim')];
% end
% ylims = [min(ylims) max(ylims)];
% xlims = [min(xlims) max(xlims)];
% 
% for i=[1 2 4 5 7 8]
%     subplot(h(i))
%     plot(cone1_xrange,cone1_uncond_rate,  '-x', 'color',[0.05 0.05 0.05], 'linewidth',1.5, 'markersize', 3)
%     axis([xlims(1)*1.05  xlims(end)*1.05 ylims(1)*0.95 ylims(2)*1.05])
%     set(gca, 'xtick', 0, 'xticklabel', '')
%     
%     if i==2 || i==5 || i == 8
%         tit1 = 'max contrast; ';
%         tmpy = [round((ylims(1)+diff(ylims)/3)*100)/100 round((ylims(1)+diff(ylims)/3*2)*100)/100];
%         tmpy = unique(tmpy);
%         set(gca, 'ytick',tmpy , 'fontsize', 6, 'fontweight', 'normal')
%     else
%         tit1 = 'min contrast; ';
%         set(gca, 'ytick', sort([0 round(ylims(2)*100)/100]), 'yticklabel', '')
%     end
%     if i<3
%         tit='y shifted';
%     elseif i<6
%         tit='x shifted';
%     else
%         tit='x AND y shifted';
%     end
%     title([tit1 tit], 'fontsize', 8, 'fontweight', 'normal')
% end
% 
% 
% subplot(h(3))
% clear t
% for i=1:size(params,1)
%     if ~isempty(params{i,1})
%         t(i,1) = params{i,1}(end);
%         t(i,2) = params1{i,1}(end);
%         t(i,3) = params2{i,1}(end);
%     end
% end
% colored_cones = zeros([size(comb), 3]);
% for i = 1:length(center_cones)
%     [a, b] = find(map==center_cones(i));
%     if i~=cone1
%         c = t(i,:)/max(t(i,:));
%     else c= [1 1 1];
%     end
%     for j = 1:length(a)
%         colored_cones(a(j),b(j), :) = c;
%     end
% end
% imagesc(colored_cones);
% axis([x y])
% set(gca, 'xtick', [],  'xticklabel', '')
% set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
% axis ij
% set(gca, 'dataaspectratio', [1 1 1])
% title({'min contrast; red: x shift better','green: y shift better, white: ref cone'}, 'fontsize', 8, 'fontweight', 'normal')
% 
% subplot(h(6))
% clear t
% for i=1:size(params,1)
%     if ~isempty(params{i,1})
%         t(i,1) = params{i,end}(end);
%         t(i,2) = params1{i,end}(end);
%         t(i,3) = params2{i,end}(end);
%     end
% end
% colored_cones = zeros([size(comb), 3]);
% for i = 1:length(center_cones)
%     [a, b] = find(map==center_cones(i));
%     if i~=cone1
%         c = t(i,:)/max(t(i,:));
%     else c= [1 1 1];
%     end
%     for j = 1:length(a)
%         colored_cones(a(j),b(j), :) = c;
%     end
% end
% imagesc(colored_cones);
% axis([x y])
% set(gca, 'xtick', [],  'xticklabel', '')
% set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
% axis ij
% set(gca, 'dataaspectratio', [1 1 1])
% title({'max contrast; red: x shift better','green: y shift better, white: ref cone'}, 'fontsize', 8, 'fontweight', 'normal')
% 
% subplot(h(9))
% colormap gray
% imagesc(sta1)
% hold on
% cnt = 1;cnt1 = 1;
% for i=1:max(map(:))
%     if ~isempty(voronoi_contours{i, 1})
%         if ~isempty(find(center_cones==i,1))
%             wid = 1.2; col = colors(cnt1, :);
%             cnt1=cnt1+1;
%         else
%             wid = 0.1;
%             col = colors2(cnt,:);
%             if cnt==5;
%                 cnt = 0;
%             end
%             cnt = cnt+1;
%         end
%         plot(voronoi_contours{i, 1}, voronoi_contours{i, 2},'color', col, 'linewidth', wid)
%     end
% end
% hold on
% if pol>0
%     cols = [0 0 0];
% else
%     cols = [1 1 1];
% end
% for i=center_cones
%     text(cones(i,2)-1,cones(i,1), int2str(i), 'color', cols, 'fontsize', 5)
% end
% axis([x y])
% set(gca, 'xtick', round([(x(1)+diff(x)/3) (x(1)+diff(x)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
% set(gca, 'ytick', round([(y(1)+diff(y)/3) (y(1)+diff(y)/3*2)]),  'fontsize', 6, 'fontweight', 'normal')
% set(gca, 'dataaspectratio', [1 1 1])
% axis ij
% title('single cone run 1',  'fontsize', 8, 'fontweight', 'normal')
% 
% 
% drawnow
% saveas(gcf,[path2save,'/tiff/', 'ID_',int2str(cellID), '_refcone_' int2str(center_cones(cone1)),'.tiff'])
% saveas(gcf,[path2save,'/svg/', 'ID_',int2str(cellID), '_refcone_' int2str(center_cones(cone1)),'.svg'])
% close(gcf)
% 
% 

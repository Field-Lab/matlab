function graphPooled(meanOnP, meanOffP, meanOnM, meanOffM,stdOnP, stdOffP, stdOnM, stdOffM, meanOnPOnM, stdOnPOnM, delta, tag)
axes_lim = [0 0.15];
hFig = figure;


plot(1.95*ones(size(stdOnP(delta(:,1)==192))), stdOnM(delta(:,1)==192), 'ko', 'MarkerFaceColor','g')
hold on
plot(ones(size(stdOnP(delta(:,1)==96))), stdOnP(delta(:,1)==96), 'ko', 'MarkerFaceColor','r')

plot(2.05*ones(size(stdOnP(delta(:,1)==48))), stdOnM(delta(:,1)==48), 'ko', 'MarkerFaceColor','b')

plot(0.95*ones(size(stdOnP(delta(:,1)==192))), stdOnP(delta(:,1)==192), 'ko', 'MarkerFaceColor','g')

plot(2.95*ones(size(stdOnPOnM(delta(:,2)==192))), stdOnPOnM(delta(:,2)==192), 'ko', 'MarkerFaceColor','g')

plot(2*ones(size(stdOnP(delta(:,1)==96))), stdOnM(delta(:,1)==96), 'ko', 'MarkerFaceColor','r')
plot(1.05*ones(size(stdOnP(delta(:,1)==48))), stdOnP(delta(:,1)==48), 'ko', 'MarkerFaceColor','b')

plot(2.95*ones(size(stdOnPOnM(delta(:,2)==192))), stdOnPOnM(delta(:,2)==192), 'ko', 'MarkerFaceColor','g')

plot(3.05*ones(size(stdOnPOnM(delta(:,2)==48))), stdOnPOnM(delta(:,2)==48), 'ko', 'MarkerFaceColor','b')
plot(3*ones(size(stdOnPOnM(delta(:,2)==96))), stdOnPOnM(delta(:,2)==96), 'ko', 'MarkerFaceColor','r')

xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'ON Midget', 'Pooled'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Variability in Speed Estimate in Different Cell Types';tag{2}});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);

s = sprintf('%c', char(176));
legend(['Speed of 42.6' s '/sec'],['Speed of 21.3' s '/sec'], ['Speed of 10.7' s '/sec'], 'location', 'northwest')


% ON Parasol std improvement by including On Midgets

axes_lim = [0 0.08];
hFig = figure;
hold on
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])

line   =zeros(2,2,length(stdOffP));
for i = 1:length(stdOnP)
    line(:,:,i) = [1, stdOnP(i); 2, stdOnPOnM(i)];
    if delta(i,2) == 48
        colors(i) = 'b';
    elseif delta(i,2) == 96
        colors(i) = 'r';
    else
        colors(i) = 'g';
    end
    
end
% colors = flipud(distinguishable_colors(length(stdOffP)));
for k = 1:length(stdOnP)
    toPlot = squeeze(line(:,:,k));
    h(k) = plot(toPlot(:,1), toPlot(:,2), 'o-','Linewidth', 2,'color',colors(k), 'MarkerFaceColor', colors(k))
end
s = sprintf('%c', char(176));
legend(h([9,3,1]),['Speed of 42.6' s '/sec'],['Speed of 21.3' s '/sec'], ['Speed of 10.7' s '/sec'], 'location', 'northeast')

xlim([0.5 2.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2]);
set(gca, 'xticklabel',{'ON Parasol', 'Pooled with ON Midgets'})
ylabel('STD/Speed')
set(gcf,'color','w');
H=sigstar({[1,2]}, 0.0452);
title({'Benefit of Pooling ON Cell Types';tag{2};'paired T test: p = 0.0452'});

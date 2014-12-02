%graph within run cell type variabilityfunc
function graphWithinRunVar(stdOnP, stdOffP, stdOnM, tag)
axes_lim = [0 0.2];
hFig = figure;
hold on
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])

line   =zeros(3,2,length(stdOffP));
for i = 1:length(stdOffP)
    line(:,:,i) = [1, stdOnP(i); 2, stdOffP(i); 3, stdOnM(i)];
end
colors = flipud(distinguishable_colors(length(stdOffP)));
for k = 1:length(stdOffP)
    toPlot = squeeze(line(:,:,k));
    plot(toPlot(:,1), toPlot(:,2), 'o-','Linewidth', 2,'color',colors(k,:), 'MarkerFaceColor', colors(k,:))
end

xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})
ylabel('STD/Speed')
set(gcf,'color','w');
title({'Variability in Cell Types within One Run';tag});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);

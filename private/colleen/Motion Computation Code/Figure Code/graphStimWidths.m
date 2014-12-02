function graphStimWidths(stdOnP, stdOffP, stdOnM, filter, tag)
%Graph stim Widths

axes_lim = [0 0.2];
hFig = figure;
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])


plot(2.95*ones(size(stdOnP(filter==8))), stdOnM(filter==8), 'ko', 'MarkerFaceColor','g')
hold on
plot(1.05*ones(size(stdOnP(filter==16))), stdOnP(filter==16), 'ko', 'MarkerFaceColor','r')

plot(0.95*ones(size(stdOnP(filter==8))), stdOnP(filter==8), 'ko', 'MarkerFaceColor','g')
plot(1.95*ones(size(stdOnP(filter==8))), stdOffP(filter==8), 'ko', 'MarkerFaceColor','g')


hold on
plot(2.05*ones(size(stdOnP(filter==16))), stdOffP(filter==16), 'ko', 'MarkerFaceColor','r')
plot(3.05*ones(size(stdOnP(filter==16))), stdOnM(filter==16), 'ko', 'MarkerFaceColor','r')


xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Effect of Stimulus Filter Width on Speed Estimate';tag});
legend('Stimulus Filter Width = 8', 'Stimulus Filter Width  = 16', 'location', 'northwest')
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);





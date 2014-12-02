function graphByDate(stdOnP, stdOffP, stdOnM, date, tag)
%Graph by experiment date





axes_lim = [0 0.2];
hFig = figure;
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])

plot(2.95*ones(size(stdOnP(date == 1))), stdOnM(date == 1), 'ko', 'MarkerFaceColor','g')
hold on
plot(1.05*ones(size(stdOnP(date == 0))), stdOnP(date == 0), 'ko', 'MarkerFaceColor','r')

plot(0.95*ones(size(stdOnP(date == 1))), stdOnP(date == 1), 'ko', 'MarkerFaceColor','g')
plot(1.95*ones(size(stdOnP(date == 1))), stdOffP(date == 1), 'ko', 'MarkerFaceColor','g')


hold on
plot(2.05*ones(size(stdOnP(date == 0))), stdOffP(date == 0), 'ko', 'MarkerFaceColor','r')
plot(3.05*ones(size(stdOnP(date == 0))), stdOnM(date == 0), 'ko', 'MarkerFaceColor','r')


xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Speed Estimate Deviation of Different Retinas';tag});
legend('2007-03-27-1', '2007-08-24-4', 'location', 'northwest')
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);





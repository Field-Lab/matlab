function graphSpeed(stdOnP, stdOffP, stdOnM, delta, tag)
axes_lim = [0 0.2];
hFig = figure;
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])
% 
% subplot_c(2,3,4)
% plot(stdOnP(delta==192), stdOffP(delta==192), 'go', 'MarkerFaceColor','g')
% 
% hold on
% plot(stdOnP(delta==96), stdOffP(delta==96), 'ro', 'MarkerFaceColor','r')
% plot(stdOnP(delta==48), stdOffP(delta==48), 'bo', 'MarkerFaceColor','b')
% 
% axis equal
% axis([axes_lim, axes_lim])
% hold on
% plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
% xlabel('ON Parasol STD/Speed')
% ylabel('OFF Parasol STD/Speed')
% 
% subplot_c(2,3,5)
% plot(stdOnP(delta==192), stdOnM(delta==192), 'go', 'MarkerFaceColor','g')
% 
% hold on
% plot(stdOnP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')
% plot(stdOnP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')
% 
% axis equal
% axis([axes_lim, axes_lim])
% hold on
% plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
% xlabel('ON Parasol STD/Speed')
% ylabel('ON Midget STD/Speed')
% 
% 
% subplot_c(2,3,6)
% plot(stdOffP(delta==192), stdOnM(delta==192), 'go', 'MarkerFaceColor','g')
% hold on
% hold on
% plot(stdOffP(delta==96), stdOnM(delta==96), 'ro', 'MarkerFaceColor','r')
% plot(stdOffP(delta==48), stdOnM(delta==48), 'bo', 'MarkerFaceColor','b')
% 
% axis equal
% axis([axes_lim, axes_lim])
% 
% hold on
% plot(axes_lim,axes_lim, '--k', 'LineWidth',2)
% xlabel('OFF Parasol STD/Speed')
% ylabel('ON Midget STD/Speed')



% subplot_c(2,3,[1 3])
title('Low Speed Variablity')
xlabel('OFF Parasol STD/Speed')
ylabel('ON Midget STD/Speed')
plot(2.95*ones(size(stdOnP(delta==192))), stdOnM(delta==192), 'ko', 'MarkerFaceColor','g')
hold on
plot(2*ones(size(stdOnP(delta==96))), stdOffP(delta==96), 'ko', 'MarkerFaceColor','r')
plot(3.05*ones(size(stdOnP(delta==48))), stdOnM(delta==48), 'ko', 'MarkerFaceColor','b')

plot(0.95*ones(size(stdOnP(delta==192))), stdOnP(delta==192), 'ko', 'MarkerFaceColor','g')
plot(1.95*ones(size(stdOnP(delta==192))), stdOffP(delta==192), 'ko', 'MarkerFaceColor','g')

plot(ones(size(stdOnP(delta==96))), stdOnP(delta==96), 'ko', 'MarkerFaceColor','r')

hold on
plot(3*ones(size(stdOnP(delta==96))), stdOnM(delta==96), 'ko', 'MarkerFaceColor','r')


hold on
plot(2.05*ones(size(stdOnP(delta==48))), stdOffP(delta==48), 'ko', 'MarkerFaceColor','b')

plot(1.05*ones(size(stdOnP(delta==48))), stdOnP(delta==48), 'ko', 'MarkerFaceColor','b')


xlim([0.5 3.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Variability in Speed Estimate in Different Cell Types';tag});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);

s = sprintf('%c', char(176));
legend(['Speed of 42.6' s '/sec'],['Speed of 21.3' s '/sec'], ['Speed of 10.7' s '/sec'], 'location', 'northwest')



function graphSpeed(meanOnP, meanOffP, meanOnM, meanOffM,stdOnP, stdOffP, stdOnM, stdOffM, delta, tag)
axes_lim = [0 0.55];
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

toKeep = stdOffM ~=0;
stdOffM = stdOffM(toKeep);
plot(4*ones(size(stdOffM(delta(toKeep)==96))), stdOffM(delta(toKeep)==96), 'ko', 'MarkerFaceColor','r')
hold on
plot(4.05*ones(size(stdOffM(delta(toKeep)==48))), stdOffM(delta(toKeep)==48), 'ko', 'MarkerFaceColor','b')
plot(3.95*ones(size(stdOffM(delta(toKeep)==192))), stdOffM(delta(toKeep)==192), 'ko', 'MarkerFaceColor','g')


xlim([0.5 4.5])
ylim(axes_lim)
set(gca, 'xtick', [1,2,3, 4]);
set(gca, 'xticklabel',{'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget'})

ylabel('STD/Speed')
set(gcf,'color','w');
title({'Variability in Speed Estimate in Different Cell Types';tag});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);

s = sprintf('%c', char(176));
legend(['Speed of 42.6' s '/sec'],['Speed of 21.3' s '/sec'], ['Speed of 10.7' s '/sec'], 'location', 'northwest')


%% Mean

figure
 set(gcf, 'Position', [280 216 1238 500]);
subplot(1,3,1)
plot(meanOnP(delta == 48), ones(size(meanOnP(delta == 48))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 48), 1.1*ones(size(meanOffP(delta == 48))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 48), 1.2*ones(size(meanOnM(delta == 48))), 'ko', 'MarkerFaceColor', 'g')
plot(10.7,0.9, 'ko', 'MarkerFaceColor', 'k')

set(gca, 'ytick' ,[0.9,1,1.1,1.2]);
set(gca, 'yticklabel', {'True Speed','ON Parasol', 'OFF Parasol', 'ON Midget'})
ylim([0.85 1.25])
xlim([10.4, 11.1])

subplot(1,3,2)
plot(meanOnP(delta == 96), ones(size(meanOnP(delta == 96))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 96), 1.1*ones(size(meanOffP(delta == 96))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 96), 1.2*ones(size(meanOnM(delta == 96))), 'ko', 'MarkerFaceColor', 'g')
plot(21.3,0.9, 'ko', 'MarkerFaceColor', 'k')
ylim([0.85 1.25])
set(gca, 'ytick' ,[]);
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
xlim([20 22.5])

subplot(1,3,3)
plot(meanOnP(delta == 192), ones(size(meanOnP(delta == 192))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 192), 1.1*ones(size(meanOffP(delta == 192))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 192), 1.2*ones(size(meanOnM(delta == 192))), 'ko', 'MarkerFaceColor', 'g')
plot(42.6,0.9, 'ko', 'MarkerFaceColor', 'k')
set(gca, 'ytick' ,[]);
xlim([39.1, 46.1])

ylim([0.85 1.25])
% xlim([40, 220])
set(gcf,'color','w');
suptitle('Speed Estimates Of Different Cell Types');

%% Equal Scale
figure
 set(gcf, 'Position', [280 216 1238 500]);
subplot(1,3,1)
plot(meanOnP(delta == 48), ones(size(meanOnP(delta == 48))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 48), 1.1*ones(size(meanOffP(delta == 48))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 48), 1.2*ones(size(meanOnM(delta == 48))), 'ko', 'MarkerFaceColor', 'g')
plot(10.7,0.9, 'ko', 'MarkerFaceColor', 'k')

set(gca, 'ytick' ,[0.9,1,1.1,1.2]);
set(gca, 'yticklabel', {'True Speed','ON Parasol', 'OFF Parasol', 'ON Midget'})
ylim([0.85 1.25])
xlim([7.2, 14.2])

subplot(1,3,2)
plot(meanOnP(delta == 96), ones(size(meanOnP(delta == 96))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 96), 1.1*ones(size(meanOffP(delta == 96))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 96), 1.2*ones(size(meanOnM(delta == 96))), 'ko', 'MarkerFaceColor', 'g')
plot(21.3,0.9, 'ko', 'MarkerFaceColor', 'k')
ylim([0.85 1.25])
set(gca, 'ytick' ,[]);
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
xlim([17.8, 24.8])

subplot(1,3,3)
plot(meanOnP(delta == 192), ones(size(meanOnP(delta == 192))), 'ko', 'MarkerFaceColor', 'r')
hold on
plot(meanOffP(delta == 192), 1.1*ones(size(meanOffP(delta == 192))), 'ko', 'MarkerFaceColor', 'b')
plot(meanOnM(delta == 192), 1.2*ones(size(meanOnM(delta == 192))), 'ko', 'MarkerFaceColor', 'g')
plot(42.6,0.9, 'ko', 'MarkerFaceColor', 'k')
set(gca, 'ytick' ,[]);
xlim([39.1, 46.1])

ylim([0.85 1.25])
% xlim([40, 220])
set(gcf,'color','w');
suptitle('Speed Estimates Shown with Equal X Scaling');






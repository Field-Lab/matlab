function graphTau(stdOnP, stdOffP, stdOnM, stdOffM, tau, tag)
%Graph by experiment date





axes_lim = [0 0.5];
hFig = figure;
% set(hFig, 'Position', [680 680 600*1.5 435*1.5])
meanOnP0_1 = mean(stdOnP(tau == 0.1));
meanOffP0_1 = mean(stdOffP(tau == 0.1));
meanOnM0_1 = mean(stdOnM(tau == 0.1));

meanOnP0_05 = mean(stdOnP(tau == 0.05));
meanOffP0_05 = mean(stdOffP(tau == 0.05));
meanOnM0_05 = mean(stdOnM(tau == 0.05));


meanOnP0_01 = mean(stdOnP(tau == 0.01));
meanOffP0_01 = mean(stdOffP(tau == 0.01));
meanOnM0_01 = mean(stdOnM(tau == 0.01));

meanOnP0_005 = mean(stdOnP(tau == 0.005));
meanOffP0_005 = mean(stdOffP(tau == 0.005));
meanOnM0_005 = mean(stdOnM(tau == 0.005));

meanOnP0_001 = mean(stdOnP(tau == 0.001));
meanOffP0_001 = mean(stdOffP(tau == 0.001));
meanOnM0_001 = mean(stdOnM(tau == 0.001));

mean0_1 = mean([meanOnP0_1, meanOffP0_1, meanOnM0_1])
mean0_05 = mean([meanOnP0_05, meanOffP0_05, meanOnM0_05])
mean0_01 = mean([meanOnP0_01, meanOffP0_01, meanOnM0_01])
mean0_005 = mean([meanOnP0_005, meanOffP0_005, meanOnM0_005])
mean0_001 = mean([meanOnP0_001, meanOffP0_001, meanOnM0_001])
% 
% plot(0.95*ones(size(stdOnP(tau == 0.1))), stdOnP(tau == 0.1), 'ko', 'MarkerFaceColor','b')
% hold on
% plot(1*ones(size(stdOnP(tau == 0.01))), stdOnP(tau == 0.01), 'ko', 'MarkerFaceColor','g')
% plot(1.05*ones(size(stdOnP(tau == 0.001))), stdOnP(tau == 0.001), 'ko', 'MarkerFaceColor','r')
% 
% plot(1.95*ones(size(stdOnP(tau == 0.1))), stdOffP(tau == 0.1), 'ko', 'MarkerFaceColor','b')
% plot(2.0*ones(size(stdOnP(tau == 0.01))), stdOffP(tau == 0.01), 'ko', 'MarkerFaceColor','g')
% plot(2.05*ones(size(stdOnP(tau == 0.001))), stdOffP(tau == 0.001), 'ko', 'MarkerFaceColor','r')
% plot(2.95*ones(size(stdOnP(tau == 0.1))), stdOnM(tau == 0.1), 'ko', 'MarkerFaceColor','b')
% plot(3*ones(size(stdOnP(tau == 0.01))), stdOnM(tau == 0.01), 'ko', 'MarkerFaceColor','g')
% plot(3.05*ones(size(stdOnP(tau == 0.001))), stdOnM(tau == 0.001), 'ko', 'MarkerFaceColor','r')
semilogx(0.1,mean0_1, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on
semilogx(0.05,mean0_05, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
semilogx(0.01, mean0_01, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
semilogx(0.005, mean0_005, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
semilogx(0.001, mean0_001, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');

xlim([0.0005 0.2])
% ylim(axes_lim)
set(gca, 'xtick', [0.001,0.005,0.01,0.05, 0.1]);
set(gca, 'xticklabel',{'0.001', '0.005', '0.01','0.05' ,'0.1'})
xlabel('\tau of Gaussian Filter')
ylabel('Average Normalized STD')
set(gcf,'color','w');
title({'Effect of Filter Width';tag});
% groups={[1.5,3]};
% H=sigstar(groups,[0.009]);

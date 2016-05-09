
% time courses
figure; hold on
default_colors = get(gca,'ColorOrder');
time = -240:8:-8;
plot(time, timecourse{1,1}, 'Color', default_colors(1,:), 'LineWidth', 2)
plot(time, timecourse{1,2}, '--', 'Color', default_colors(1,:))
time = 2*(-240:8:-8);
plot(time, timecourse{2,1}, 'Color', default_colors(2,:), 'LineWidth', 2)
plot(time, timecourse{2,2}, '--', 'Color', default_colors(2,:))
plot(time, timecourse{3,1}, 'Color', default_colors(5,:), 'LineWidth', 2)
plot(time, timecourse{3,2}, '--', 'Color', default_colors(5,:))
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')

% RGB balance
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(RGB{1,1}, 'Color', default_colors(1,:), 'LineWidth', 2)
plot(RGB{1,2}, '--', 'Color', default_colors(1,:))
plot(RGB{2,1},'Color', default_colors(2,:), 'LineWidth', 2)
plot(RGB{2,2}, '--', 'Color', default_colors(2,:))
plot(RGB{3,1},'Color', default_colors(5,:), 'LineWidth', 2)
plot(RGB{3,2}, '--', 'Color', default_colors(5,:))
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
xlim([0 4])
set(gca,'XTick',1:3)
set(gca,'XTickLabels',{'R', 'G', 'B'})
ylim([0 1])
title('RGB Balance')

% RF size
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(2*ones(length(sta_mean{1,1}),1),sta_mean{1,1}*8,'.','Color', default_colors(1,:))
plot(2.5*ones(length(sta_mean{1,2}),1),sta_mean{1,2}*8,'.','Color', default_colors(1,:))
plot(3*ones(length(sta_mean{2,1}),1),sta_mean{2,1}*10,'.','Color', default_colors(2,:))
plot(3.5*ones(length(sta_mean{2,2}),1),sta_mean{2,2}*10,'.','Color', default_colors(2,:))
plot(4*ones(length(sta_mean{3,1}),1),sta_mean{3,1}*10,'.','Color', default_colors(5,:))
plot(4.5*ones(length(sta_mean{3,2}),1),sta_mean{3,2}*10,'.','Color', default_colors(5,:))
xlim([1.5 5])
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('RF size in screen pixels')

% spike rate
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(2*ones(length(spike_rate{1,1}),1),spike_rate{1,1},'.','Color', default_colors(1,:))
plot(2.5*ones(length(spike_rate{1,2}),1),spike_rate{1,2},'.','Color', default_colors(1,:))
plot(3*ones(length(spike_rate{2,1}),1),spike_rate{2,1},'.','Color', default_colors(2,:))
plot(3.5*ones(length(spike_rate{2,2}),1),spike_rate{2,2},'.','Color', default_colors(2,:))
plot(4*ones(length(spike_rate{3,1}),1),spike_rate{3,1},'.','Color', default_colors(5,:))
plot(4.5*ones(length(spike_rate{3,2}),1),spike_rate{3,2},'.','Color', default_colors(5,:))
xlim([1.5 5])
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('spike rate')

clear time

%%
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(2*ones(length(bps_exp{1,1}),1),bps_exp{1,1},'.','Color', default_colors(1,:))
plot(2.5*ones(length(bps_exp{1,2}),1),bps_exp{1,2},'.','Color', default_colors(1,:))
plot(3*ones(length(bps_exp{2,1}),1),bps_exp{2,1},'.','Color', default_colors(2,:))
plot(3.5*ones(length(bps_exp{2,2}),1),bps_exp{2,2},'.','Color', default_colors(2,:))
plot(4*ones(length(bps_exp{3,1}),1),bps_exp{3,1},'.','Color', default_colors(5,:))
plot(4.5*ones(length(bps_exp{3,2}),1),bps_exp{3,2},'.','Color', default_colors(5,:))
xlim([1.5 5])
ylim([-0.05, 0.8])
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('BPS')

% mu
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(2*ones(length(mu_exp{1,1}),1),mu_exp{1,1},'.','Color', default_colors(1,:))
plot(2.5*ones(length(mu_exp{1,2}),1),mu_exp{1,2},'.','Color', default_colors(1,:))
plot(3*ones(length(mu_exp{2,1}),1),mu_exp{2,1},'.','Color', default_colors(2,:))
plot(3.5*ones(length(mu_exp{2,2}),1),mu_exp{2,2},'.','Color', default_colors(2,:))
plot(4*ones(length(mu_exp{3,1}),1),mu_exp{3,1},'.','Color', default_colors(5,:))
plot(4.5*ones(length(mu_exp{3,2}),1),mu_exp{3,2},'.','Color', default_colors(5,:))
xlim([1.5 5])
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('mu')

% time courses
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(K_time_exp{1,1}/max(K_time_exp{1,1}), 'Color', default_colors(1,:), 'LineWidth', 2)
plot(K_time_exp{1,2}/max(K_time_exp{1,2}), '--', 'Color', default_colors(1,:))
plot(K_time_exp{2,1}/max(K_time_exp{2,1}), 'Color', default_colors(2,:), 'LineWidth', 2)
plot(K_time_exp{2,2}/max(K_time_exp{2,2}), '--', 'Color', default_colors(2,:))
plot(K_time_exp{3,1}/max(K_time_exp{3,1}), 'Color', default_colors(5,:), 'LineWidth', 2)
plot(K_time_exp{3,2}/max(K_time_exp{3,2}), '--', 'Color', default_colors(5,:))
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('GLM Temporal filter')

% PS
figure; hold on
default_colors = get(gca,'ColorOrder');
plot(PS_exp{1,1}, 'Color', default_colors(1,:), 'LineWidth', 2)
plot(PS_exp{1,2}, '--', 'Color', default_colors(1,:))
plot(PS_exp{2,1}, 'Color', default_colors(2,:), 'LineWidth', 2)
plot(PS_exp{2,2}, '--', 'Color', default_colors(2,:))
plot(PS_exp{3,1}, 'Color', default_colors(5,:), 'LineWidth', 2)
plot(PS_exp{3,2}, '--', 'Color', default_colors(5,:))
legend('2012-08-09-3 ON','2012-08-09-3 OFF', '2012-09-27-3 ON','2012-09-27-3 OFF', '2013-08-19-6 ON', '2013-08-19-6 OFF')
title('GLM PS filter')
toc

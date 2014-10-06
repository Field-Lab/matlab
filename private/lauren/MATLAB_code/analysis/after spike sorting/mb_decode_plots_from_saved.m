clear all

visWS = load('/snle/lab/Experiments/Array/Analysis/2011-06-24-5/vis_mb_decode_workspace.mat');
elecWS = load('/snle/lab/Experiments/Array/Analysis/2011-06-24-5/elec_mb_decode_workspace.mat');

trial_to_highlight = 70;

bar_speed = visWS.bar_speed;

%% motion curve, single vis trail

figure; hold on
plot(bar_speed*[1 1], [0 max(visWS.sig_strength_curve{trial_to_highlight})*1.2], 'k--')
plot(visWS.run_opt.vel, visWS.sig_strength_curve{trial_to_highlight})
xlabel('speed (pixels/s)')
ylabel('signal strength')
title('speed tuning, single visual response trial')

%% plot motion signal curves (vis) on top of each other with chosen trial highlighted

%lineColors = lines(nTrials);
figure; hold on
plot(bar_speed*[1 1], [0 max(cellfun(@max, visWS.sig_strength_curve))]*1.1, 'k--')
for jj = 1:visWS.nTrials
    plot(visWS.run_opt.vel, visWS.sig_strength_curve{jj}, 'color', [0.8 0.8 0.8])
end
plot(visWS.run_opt.vel, visWS.sig_strength_curve{trial_to_highlight}, 'r-', 'linewidth', 1)
plot(visWS.vel_ests_brute(trial_to_highlight), max(visWS.sig_strength_curve{trial_to_highlight}), 'r.');

xlabel('speed (pixels/second)')
ylabel('net motion signal')
title('speed tuning curve of each trial (visual stimulation)')


%% plot all motion curves, electrical stimulation trials with target visual response

figure; hold on
plot(bar_speed*[1 1], [0 max(cellfun(@max, visWS.sig_strength_curve))*1.1], 'k--')

for jj = 1:elecWS.nTrials
    plot(elecWS.run_opt.vel, elecWS.sig_strength_curve{jj}, 'color', [0.8 0.8 0.8])
end

%target visual response
plot(visWS.run_opt.vel, visWS.sig_strength_curve{trial_to_highlight}, 'r-', 'linewidth', 1)
plot(visWS.vel_ests_brute(trial_to_highlight), max(visWS.sig_strength_curve{trial_to_highlight}), 'r.');

xlabel('speed (pixels/second)')
ylabel('net motion signal')
title('speed tuning curve of each trial (electrical stimulation and target visual response)')



%% plot all motion curves, visual and electrical stimulation trials

figure('position', [100 200 400 600]);
yLimCurve = [-0.2 12];
yLimHist = [0 12];

axes('position', [0.1 0.55 0.8 0.35]); hold on
%plot(bar_speed*[1 1], [0 max(cellfun(@max, visWS.sig_strength_curve))*1.1], 'k--')
plot(bar_speed*[1 1], yLimCurve, 'k--')

for jj = 1:visWS.nTrials
    plot(visWS.run_opt.vel, visWS.sig_strength_curve{jj}, 'color', [0.5 0.5 0.5])
end
for jj = 1:elecWS.nTrials
    plot(elecWS.run_opt.vel, elecWS.sig_strength_curve{jj}, 'color', [0.8 0 0])
end

plot(visWS.run_opt.vel, visWS.sig_strength_curve{trial_to_highlight}, 'k-', 'linewidth', 1)
plot(visWS.vel_ests_brute(trial_to_highlight), max(visWS.sig_strength_curve{trial_to_highlight}), 'k.');
set(gca, 'ylim', yLimCurve, 'xlim', [visWS.run_opt.vel(1) visWS.run_opt.vel(end)])

xlabel('speed (pixels/second)')
ylabel('net motion signal')
title(['speed tuning curve of each trial:' 10 'electrical stim (red), visual stim (gray) and target visual response (black)'])

axes('position', [0.1 0.3 0.8 0.18]); hold on
plot(bar_speed*[1 1], yLimHist, 'k--')
histCenters = visWS.run_opt.vel;
hist(visWS.vel_ests_brute, histCenters);
h = findobj(gca, 'type', 'patch');
set(h, 'facecolor', [0.5 0.5 0.5], 'edgecolor', 'none');
set(gca, 'ylim', yLimHist, 'xlim', [visWS.run_opt.vel(1) visWS.run_opt.vel(end)])

axes('position', [0.1 0.05 0.8 0.18]); hold on
plot(bar_speed*[1 1], yLimHist, 'k--')
hist(elecWS.vel_ests_brute, histCenters);
h = findobj(gca, 'type', 'patch');
%h = h(cellfun(@(x) strcmpi('flat', x), get(h, 'facecolor')));
set(h, 'facecolor', [0.8 0 0], 'edgecolor', 'none');
set(gca, 'ylim', yLimHist, 'xlim', [visWS.run_opt.vel(1) visWS.run_opt.vel(end)])

xlabel('speed estimate')
ylabel('number of trials')

%% plot histogram' stim_type '

figure
histCenters = 200:300;
hist(vel_ests_brute, histCenters);
title(['speed estimates (brute force, ' stim_type ' stimulation)'])


%% plot histogram

figure
histCenters = 200:300;
hist(vel_ests_eff, histCenters);
title(['speed estimates (search-based, ' stim_type ' stimulation)'])
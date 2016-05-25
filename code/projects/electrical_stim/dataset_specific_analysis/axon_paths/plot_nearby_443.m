figure;

scatter(elec_dist(1:5,2),elec_dist(1:5,3), 100, '*'); % Too noisy
hold on;
over_thresh = [0 5;120 5];
plot(over_thresh(:,1),over_thresh(:,2), '--');
hold on;
scatter(elec_dist(6:10,2),elec_dist(6:10,3), 100,'filled'); % Just right
hold on;
scatter(elec_dist(11:13,2),elec_dist(11:13,3), 100, 'filled'); % Over Threshold


legend('Too noisy to analyze', 'Cell not activated by stim electrode');


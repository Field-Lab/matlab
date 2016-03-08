figure;
hold on;

scatter(Amps(Activity == 0,1), Amps(Activity == 0,2), 'b', 'filled');
scatter(Amps(Activity == 1,1), Amps(Activity == 1,2), 'r', 'filled');

xlabel('Electrode 17 (uA)');
ylabel('Electrode 372 (uA)');

% plot dashed lines

plot(0:0.1:3.5, 1.31.*ones(size(0:0.1:3.5)), 'k--');
plot(1.41*ones(size(0:0.1:3.5)), 0:0.1:3.5, 'k--');
clear
%% datarun 1 = class, 2 = full_rep, 3 = vertical half, 4 = left third, 5 = left two thirds
datarun{1} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data022-from-data018_data019_data020_data021_data022/data022-from-data018_data019_data020_data021_data022');
datarun{2} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data020-from-data018_data019_data020_data021_data022/data020-from-data018_data019_data020_data021_data022');
datarun{3} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data021-from-data018_data019_data020_data021_data022/data021-from-data018_data019_data020_data021_data022');
datarun{4} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data018-from-data018_data019_data020_data021_data022/data018-from-data018_data019_data020_data021_data022');
datarun{5} = load_data('/Volumes/Analysis/2015-08-17-2/data018-data022/data019-from-data018_data019_data020_data021_data022/data019-from-data018_data019_data020_data021_data022');

datarun{1} = load_params(datarun{1});
for i = 1:5
    datarun{i} = load_neurons(datarun{i});
end

plot_rf_fit(datarun{1}, 'Off Parasol')

%%
cid = get_cell_indices(datarun{1}, 'Off Parasol');
% MSE = zeros(length(cid), 1);
loc = zeros(length(cid), 2);
for i_cell = 10%length(cid)
    for i_run = 2:5
        spikes_concat = [];
        spikes = datarun{i_run}.spikes{cid(i_cell)};
        trial_starts = datarun{i_run}.triggers([true; diff(datarun{i_run}.triggers)>0.9]);
        figure; hold on
        for i_trial = 1:length(trial_starts)
            trial_spikes = spikes( (spikes>trial_starts(i_trial)) & (spikes<(trial_starts(i_trial) + 15)) ) - trial_starts(i_trial);
            plot(trial_spikes, i_trial*ones(size(trial_spikes)), 'k.')
            spikes_concat = [spikes_concat; trial_spikes];
        end
        hold off
        spikes_concat = sort(spikes_concat);
        PSTH_temp = zeros(15000,1);
        PSTH_temp(ceil(spikes_concat*1000)) = 1;
        PSTH{i_run} = conv(PSTH_temp, gausswin(1000), 'same');
        % plot(PSTH{i_run})
    end
    if ~mod(i_cell, 10); disp(i_cell); end
    NSEM_Corr_vert(i_cell) = err(PSTH{2}, PSTH{3});
    NSEM_Corr_half(i_cell) = err(PSTH{2}, PSTH{4});
    NSEM_Corr_twothirds(i_cell) = err(PSTH{2}, PSTH{5});
    loc(i_cell, :) = datarun{1}.vision.sta_fits{cid(i_cell)}.mean;
end

%%
figure; plot(loc(:,2), NSEM_Corr_vert, '.', 'MarkerSize', 10)
hold on; plot([20 20], [-0.2 1.3], 'LineWidth', 2)
figure; plot(loc(:,1), NSEM_Corr_half, '.', 'MarkerSize', 10)
hold on; plot([27 27], [-0.2 1.3], 'LineWidth', 2)
figure; plot(loc(:,1), NSEM_Corr_twothirds, '.', 'MarkerSize', 10)
hold on; plot([50 50], [-0.2 1.3], 'LineWidth', 2)


%%
clear
% datarun 1 = class, 2 = full_rep, 3 = left third, 4 = left two thirds
datarun{1} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data013-from-data013_data022_data023_data024/data013-from-data013_data022_data023_data024');
datarun{2} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data022-from-data013_data022_data023_data024/data022-from-data013_data022_data023_data024');
datarun{3} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data023-from-data013_data022_data023_data024/data023-from-data013_data022_data023_data024');
datarun{4} = load_data('/Volumes/Analysis/2015-08-17-6/data022-24-repeats/data024-from-data013_data022_data023_data024/data024-from-data013_data022_data023_data024');

datarun{1} = load_params(datarun{1});
for i = 1:4
    datarun{i} = load_neurons(datarun{i});
end

figure; 
plot_rf_fit(datarun{1}, 'Off Midget')

%%
cid = get_cell_indices(datarun{1}, 'Off Midget');
% MSE = zeros(length(cid), 1);
loc = zeros(length(cid), 2);
for i_cell = 1:length(cid)
    for i_run = 2:4
        spikes_concat = [];
        spikes = datarun{i_run}.spikes{cid(i_cell)};
        trial_starts = datarun{i_run}.triggers([true; diff(datarun{i_run}.triggers)>0.9]);
        for i_trial = 1:length(trial_starts)
            spikes_concat = [spikes_concat; spikes( (spikes>trial_starts(i_trial)) & (spikes<(trial_starts(i_trial) + 15)) ) - trial_starts(i_trial)];
        end
        spikes_concat = sort(spikes_concat);
        PSTH_temp = zeros(15000,1);
        PSTH_temp(ceil(spikes_concat*1000)) = 1;
        PSTH{i_run} = conv(PSTH_temp, gausswin(500), 'same');
    end
    disp(i_cell)
    % NSEM_Corr_vert(i_cell) = corr(PSTH{2}, PSTH{3});
    NSEM_Corr_half(i_cell) = err(PSTH{2}, PSTH{3});
    NSEM_Corr_twothirds(i_cell) = err(PSTH{2}, PSTH{4});
    loc(i_cell, :) = datarun{1}.vision.sta_fits{cid(i_cell)}.mean;
end

%%
% figure; plot(loc(:,2), NSEM_Corr_vert, '.')
% hold on; plot([20 20], [-0.2 1.3])
figure; plot(loc(:,1), NSEM_Corr_half, '.', 'MarkerSize', 10)
hold on; plot([40 40], [-0.2 1.3], 'LineWidth', 2)
figure; plot(loc(:,1), NSEM_Corr_twothirds, '.', 'MarkerSize', 10)
hold on; plot([50 50], [-0.2 1.3], 'LineWidth', 2)

%%
clear
% datarun 1 = class, 2 = full_rep, 3 = left third, 4 = left two thirds
datarun{1} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data017-from-data012_data013_data017/data017-from-data012_data013_data017');
datarun{2} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data012-from-data012_data013_data017/data012-from-data012_data013_data017');
datarun{3} = load_data('/Volumes/Analysis/2015-09-23-0/data012-data013-data017-norefit/data013-from-data012_data013_data017/data013-from-data012_data013_data017');

datarun{1} = load_params(datarun{1});
for i = 1:3
    datarun{i} = load_neurons(datarun{i});
end

%%
figure; 
plot_rf_fit(datarun{1}, 'Off Parasol')
hold on; plot_rf_fit(datarun{1}, 67, 'fill', true)
hold on; plot_rf_fit(datarun{1}, 4998, 'fill', true)


%%
cid = get_cell_indices(datarun{1}, 'Off Parasol');
% MSE = zeros(length(cid), 1);
clear NSEM_Corr_half
loc = zeros(length(cid), 2);
for i_cell = 2%1:length(cid)
    disp(cid(i_cell))
    for i_run = 2:3
        spikes_concat = [];
        spikes = datarun{i_run}.spikes{cid(i_cell)};
        trial_starts = datarun{i_run}.triggers([true; diff(datarun{i_run}.triggers)>0.9]);
        for i_trial = 1:length(trial_starts)
            spikes_concat = [spikes_concat; spikes( (spikes>trial_starts(i_trial)) & (spikes<(trial_starts(i_trial) + 15)) ) - trial_starts(i_trial)];
        end
        spikes_concat = sort(spikes_concat);
        PSTH_temp = zeros(15000,1);
        PSTH_temp(ceil(spikes_concat*1000)) = 1;
        PSTH{i_run} = conv(PSTH_temp, gausswin(500), 'same');
    end
    plot(PSTH{2})
    hold on; plot(PSTH{3}); hold off;
    disp(err(PSTH{2}, PSTH{3}))
     %pause();
    disp(i_cell)
    %NSEM_Corr_vert(i_cell) = corr(PSTH{2}, PSTH{3});
    NSEM_Corr_half(i_cell) = err(PSTH{2}, PSTH{3});
    loc(i_cell, :) = 5.5*8*datarun{1}.vision.sta_fits{cid(i_cell)}.mean;
end


%%
% figure; plot(loc(:,2), NSEM_Corr_vert, '.')
% hold on; plot([20 20], [-0.2 1.3])
figure; plot(loc(:,2), NSEM_Corr_half/15000, '.', 'MarkerSize', 10)
hold on; plot([5.5*8*20 5.5*8*20], [-0.2 0.7], 'LineWidth', 2)
xlim([0 2000])
ylim([-0.25 0.75])
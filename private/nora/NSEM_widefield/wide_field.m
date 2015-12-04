function [NSEM_Corr, loc, avg_profile, PSTH_ex_cells] = wide_field(dataruns, cell_type, stim_end, example_cells)

%% datarun 1 = class, 2 = full_rep, 3 = split

dataruns{1} = load_params(dataruns{1});
for i = 1:length(dataruns)
    dataruns{i} = load_neurons(dataruns{i});
end

%% Plot mosaic
% if example_cells
% figure; plot_rf_fit(dataruns{1}, cell_type)
%     ex_cells = get_cell_indices(dataruns{1}, example_cells);
%     PSTH_ex_cells = cell(length(example_cells), length(dataruns)-1);
%     for i=1:length(example_cells)
%         hold on; plot_rf_fit(dataruns{1}, example_cells(i), 'fill', true)
%     end
% end

% %% avg RF
% dataruns{1} = load_sta(dataruns{1});
% avg_rf = sum(get_average_rf(dataruns{1}, cell_type, 'scale', 5),3);
% 
% %
% n_angles = 50;
% slice = zeros(200,1);
% for i = 1:n_angles
%     angle = i*360/n_angles;
%     temp = imrotate(avg_rf, angle, 'crop');
%     slice = slice+ improfile(temp, [200 200], [1 200]);
% end
% slice = slice/n_angles;
% a = fit((1:200)'/5,slice, fittype('gauss1'));
% width = a.c1/sqrt(2);
% avg_profile = [slice, (1:200)'/(5*width)];
avg_profile = 0;
width = 1;

%% Organize PSTH and calculate errors
cid = get_cell_indices(dataruns{1}, cell_type);
loc = zeros(length(cid),1);
NSEM_Corr = zeros(length(cid), length(dataruns)-1);
i_ex_cell = 0;
for i_cell = 1:length(cid)
    PSTH = cell(1, length(dataruns)-1);
    for i_run = 2:length(dataruns)
        spikes_concat = [];
        spikes = dataruns{i_run}.spikes{cid(i_cell)};
        trial_starts = dataruns{i_run}.triggers([true; diff(dataruns{i_run}.triggers)>0.9]);
        for i_trial = 1:length(trial_starts)
            trial_spikes = spikes( (spikes>trial_starts(i_trial)) & (spikes<(trial_starts(i_trial) + 15)) ) - trial_starts(i_trial);
            spikes_concat = [spikes_concat; trial_spikes];
        end
        spikes_concat = ceil(spikes_concat*100);
        PSTH_temp = zeros(1500,1);
        for i = 1:1500
            PSTH_temp(i) = sum(spikes_concat == i);
        end
        PSTH{i_run} = 100*conv(PSTH_temp, gausswin(5), 'same')/(40*sum(gausswin(5))); % in hertz
        NSEM_Corr(i_cell,i_run) = err(PSTH{2}, PSTH{i_run})/1500;
    end
    %if any(example_cells) && any(ex_cells == cid(i_cell))
    loc(i_cell) = (dataruns{1}.vision.sta_fits{cid(i_cell)}.mean(2)-stim_end(1))/width;
    if loc < -4
        i_ex_cell = i_ex_cell + 1;
        PSTH_ex_cells{i_ex_cell} = PSTH;
    end
    if ~mod(i_cell, 50); disp(i_cell); end
end

% 
% if length(stim_end) == 2
%     loc = [loc; loc-diff(stim_end)/width];
% end



% temp
% plot(loc, [NSEM_Corr(:,3)], '.')
% text(loc+0.1, [NSEM_Corr(:,3)], [strsplit(num2str(cid))]);

end
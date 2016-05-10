function [NSEM_Corr, loc, avg_profile, PSTH_ex_cells] = wide_field_newnorm(dataruns, cell_type, stim_end, example_cells)

%% datarun 1 = class, 2 = full_rep, 3 = split
modu=[];
no_modu= [];

dataruns{1} = load_params(dataruns{1});
for i = 1:length(dataruns)
    dataruns{i} = load_neurons(dataruns{i});
end

%% Plot mosaic
if example_cells
figure; plot_rf_fit(dataruns{1}, cell_type)
    ex_cells = get_cell_indices(dataruns{1}, example_cells);
    PSTH_ex_cells = cell(length(example_cells), length(dataruns)-1);
    for i=1:length(example_cells)
        hold on; plot_rf_fit(dataruns{1}, example_cells(i), 'fill', true)
    end
end

%% avg RF
dataruns{1} = load_sta(dataruns{1});
avg_rf = sum(get_average_rf(dataruns{1}, cell_type, 'scale', 5),3);

n_angles = 10;
slice = zeros(201,1);
range = -100:100;

for i = 1:4
    figure(1); imagesc(avg_rf); axis image; title(i)
    prof = improfile;
    prof = [zeros(100,1); prof; zeros(100,1)];
    [~,center] = max(prof);
    slice = slice + prof(center+range);
end
slice = slice/n_angles;
plot(slice)

a = fit((1:201)'/5,slice, fittype('gauss1'));
width = a.c1/sqrt(2);
avg_profile = [slice, (1:201)'/(5*width)];
% avg_profile = 0;
% width = 1;

%% Organize PSTH and calculate errors
cid = get_cell_indices(dataruns{1}, cell_type);
cells = get_cell_ids(dataruns{1}, cell_type);
for i_run = 2:length(dataruns)
        prepped_data{i_run} = interleaved_data_prep(dataruns{i_run}, 1800, 40, 'cell_spec', cells, 'visual_check', 0);
end
loc = zeros(length(cid),1);
NSEM_Corr = zeros(length(cid), length(dataruns)-1);
i_ex_cell = 0;
for i_cell = 1:length(cid)
    PSTH = cell(1, length(dataruns)-1);
    for i_run = 2:length(dataruns)
        %         spikes_concat = [];
        %         spikes = dataruns{i_run}.spikes{cid(i_cell)};
        %         trial_starts = dataruns{i_run}.triggers([true; diff(dataruns{i_run}.triggers)>0.9]);
        %         for i_trial = 1:length(trial_starts)
        %             trial_spikes = spikes( (spikes>trial_starts(i_trial)) & (spikes<(trial_starts(i_trial) + 15)) ) - trial_starts(i_trial);
        %             spikes_concat = [spikes_concat; trial_spikes];
        %         end
        %         spikes_concat = ceil(spikes_concat*100);
        %         PSTH_temp = zeros(1500,1);
        %         for i = 1:1500
        %             PSTH_temp(i) = sum(spikes_concat == i);
        %         end
        %         PSTH{i_run} = 100*conv(PSTH_temp, gausswin(5), 'same')/(40*sum(gausswin(5))); % in hertz
        PSTH{i_run} = IDP_plot_PSTH(prepped_data{i_run}, i_cell, 0, 1/120, 10);
        NSEM_Corr(i_cell,i_run) = err(PSTH{2}, PSTH{i_run});
    end
    loc(i_cell) = (dataruns{1}.vision.sta_fits{cid(i_cell)}.mean(2)-stim_end(1))/width;
    
    if any(example_cells) && any(ex_cells == cid(i_cell))
        i_ex_cell = i_ex_cell + 1;
        PSTH_ex_cells{i_ex_cell} = PSTH;
    end
%    if loc(i_cell) < -3
%         
%         f = figure(1);
%         set(f, 'Position', [100 100 600 200]);
%         plot(PSTH{3})
%         hold on
%         for i = 1:15
%             plot([i i]*100, [0 100], 'k')
%         end
%         hold off
%         ylim([0 max(PSTH{3})]);
%         title(loc(i_cell))
%         pause()
% %         if ~x %if click
% %             modu = [modu loc(i_cell)];
% %             close(f)
% %             disp('mod')
% %         elseif x % if hit enter
% %             no_modu = [no_modu loc(i_cell)];
% %             close(f)
% %             disp('no mod')
% %         end
%     end
%    modu = 0;
%    no_modu = 0;
    if ~mod(i_cell, 50); disp(i_cell); end
end

% 
if length(stim_end) == 2
    loc = [loc; loc-diff(stim_end)/width];
elseif length(stim_end) == 3
    loc = loc*width+stim_end(1);
    loc = [(loc-stim_end(1))/width; (loc-stim_end(2))/width; (loc-stim_end(3))/width];
end


% temp
% plot(loc, [NSEM_Corr(:,3)], '.')
% text(loc+0.1, [NSEM_Corr(:,3)], [strsplit(num2str(cid))]);

end
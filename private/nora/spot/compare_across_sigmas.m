 clear 
params_201602171
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

n_subgroups = length(cell_idx);

%%
%%{
% load reg repeats
for i_reg = 1:n_reg
    datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
    datarun{i_reg} = load_neurons(datarun{i_reg});
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 0);
    clear datarun
end
%}

%%{
    % plot reg repeats for stability check
    for i_cell = 1:size(reg_data{1}.testspikes,2)
        % figure; hold on
        for i_reg = 1:n_reg
            reg_PSTH{i_cell} = IDP_plot_PSTH(reg_data{i_reg}, i_cell, 'color', 0);
        end
        title(['cell ' num2str(i_cell)]);
    end
%}
%close all

%%
datarun = load_data([ Analysis_Path masking '/' masking]);
datarun = load_neurons(datarun);
mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
clear datarun

% separate each condition
idx = 1:30;
for count = 1:(n_masks*2)
    condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
end

%%
% plot the PSTH for each condition for each cell
start = 1;
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        full_firing = sum(reg_PSTH{i_cell});
        for i = 1:length(mask_conditions{subgroup})
            mask = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            %comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 100);
            %for fixation = 1:23
            Firing(i_cell, i) = sum(mask)/full_firing;
                %start = interval_frame(fixation);
            %end
            %start = 1;
            %Firing(i_cell, i, 2) = sum(comp);
        end
        %figure(1); hold on;
        %plot(mask{1}, reg_PSTH{i_cell}, '.k')
    end
end



%hold on; plot([0 300], [0 300])
%        pause()
% plot([2 4 5 6], Firing(:,1))
save([fig_save '/Firing.mat'], 'Firing')

%%
%{
for i = 1:100
    load(['/Volumes/Lab/Users/Nora/new_stim_nora/NSEM intervals/matfiles/movie_chunk_' num2str(i) '.mat'])
    interval(i) = size(movie_chunk, 3);
end
default_colors = get(gca,'ColorOrder');
interval_frame = cumsum(interval);
%}

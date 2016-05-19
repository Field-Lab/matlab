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
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 1);
    clear datarun
end
%}

%{
    % plot reg repeats for stability check
    for i_cell = 1:size(reg_data{1}.testspikes,2)
        figure; hold on
        for i_reg = 1:n_reg
            IDP_plot_PSTH(reg_data{i_reg}, i_cell, 'color', i_reg);
        end
        title(['cell ' num2str(i_cell)]);
    end
%}
close all

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
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        [reg, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0);
        for i = 1:length(mask_conditions{subgroup})
            color = [1 1 1]*0;%(1-sigmas(mask_conditions{subgroup}(i))/12);
            figure(i); hold on
            mask = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0);
            comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0);
            Corr(i_cell, i, 1) = corr(reg, mask);
            Corr(i_cell, i, 2) = corr(reg, comp);
            Corr(i_cell, i, 3) = corr(comp, mask);
            %plot(time, reg, 'k');
            %exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_sigma_' num2str(sigmas(mask_conditions{subgroup}(i)))], 'Bounds', 'loose', 'Color', 'rgb')
            %clf
        end
    end
end
close all
save([fig_save '/Corr.mat'], 'Corr')